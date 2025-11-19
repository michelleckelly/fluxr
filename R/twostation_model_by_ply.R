twostation_model_by_ply <- function(model_fun, data, day_start, day_end){
  data.plys <- as.data.frame(data)
  data.plys <- arrange(data.plys, solar.time)
  
  min_timestep <- 
    streamMetabolizer::mm_get_timestep(unique(data$solar.time), 
                                       format = "unique")
  if(length(min_timestep) == 1 && !is.na(min_timestep)[1] && min_timestep <= 0){
    timesteps <- as.numeric(diff(unique(data$solar.time)), units = "days")
    timegoof <- which.min(timesteps) + c(0,1)
    stop("min timestep is <= 0: ", format(min_timestep, digits=3), " days from ",
         data$solar.time[timegoof[1]], " (time ", timegoof[1], ") to ",
         data$solar.time[timegoof[2]], " (time ", timegoof[2], ")")
  }
  if(!is.null(data_daily)){
    min_datestep <- 
      streamMetabolizer::mm_get_timestep(unique(data$date), 
                                         format = "unique")
    if(length(min_datestep) > 0 && !is.na(min_datestep)[1] && min_datestep[1] <= 0){
      timesteps <- as.numeric(diff(data_daily$date), units="days")
      timegoof <- which.min(timesteps) + c(0,1)
      stop("min datestep is <= 0: ", min_datestep, " days from ",
           data_daily$date[timegoof[1]], " (row ", timegoof[1], ") to ",
           data_daily$date[timegoof[2]], " (row ", timegoof[2], ")")
    }
  }
  doyhr <- convert_date_to_doyhr(data.plys$solar.time)
  data.plys$date <- as.character(as.Date(
    lubridate::floor_date(data.plys$solar.time, "year") + 
      as.difftime(floor(doyhr)-1, units = "days")))
  data.plys$hour <- 24*(doyhr %% 1)
  
  # Add in a 1 second fudge factor
  day_start_fudge <- day_start - 1/(60*60*24)
  day_end_fudge <- day_end - 1/(60*60*24)
  
  get_runs <- function(date.group) {
    values <- start <- end <- '.dplyr.var'
    replace(date.group, is.na(date.group), '') %>%
      rle() %>%
      unclass %>%
      as_tibble %>%
      mutate(end=cumsum(lengths), start=c(1, 1+end[-n()])) %>%
      filter(values != '') %>%
      select(date=values, start, end)
  }
  
  # streammetabolizer function labels dates as odd or even 
  date_rows <- get_runs(data.plys$date)
  odd.dates <- date_rows[which(seq_len(nrow(date_rows)) %% 2 == 1),]$date
  even.dates <- date_rows[which(seq_len(nrow(date_rows)) %% 2 == 0),]$date
  
  # Assign each data row to one or two date plys
  dt.NA <- as.character(NA)
  if(nrow(data.plys) == 0){
    data.plys <- bind_cols(
      data.plys,
      tibble::tibble(odd.date.group = dt.NA, even.date.group = dt.NA)[c(),]
    )
    out_list <- list()
  } else{
    data.plys <- bind_cols(
      data.plys,
      bind_rows(lapply(seq_len(nrow(date_rows)), function(dt){
        run <- date_rows[dt,]
        dt.today <- run$date
        dt.yesterday <- if(dt>1) date_rows[[dt-1, "date"]] else dt.NA
        dt.tomorrow <- if(dt<nrow(date_rows)) date_rows[[dt+1, "date"]] else dt.NA
        data.plys.rows <- run$start:run$end
        hr <- data.plys[data.plys.rows, "hour"]
        primary.date <- c(dt.today, dt.NA)[ifelse(hr >= day_start_fudge & hr < day_end_fudge, 1, 2)]
        secondary.date <- c(dt.yesterday, dt.NA, dt.tomorrow)[ifelse(hr <= (day_end_fudge - 24), 1, ifelse(hr < (24 + day_start_fudge), 2, 3))]
        if(dt.today %in% odd.dates){
          tibble(odd.date.group = primary.date, even.date.group = secondary.date)
        } else{
          tibble(odd.date.group = secondary.date, even.date.group = primary.date)
        }
      }))) %>%
      as.data.frame(stringsAsFactors = FALSE)
     
    date.pairings <- setNames(unique(data.plys[, c("odd.date.group", "even.date.group")]), c("a", "b"))
    date.pairings <- rbind(date.pairings, setNames(date.pairings, c("b", "a")))
    date.pairings <- date.pairings[!is.na(date.pairings$a),]
    solo.dates <- date.pairings[is.na(date.pairings$b),]$a
    tbl.dates <- table(date.pairings$a)
    double.dates <- names(tbl.dates)[tbl.dates > 1]
    unique.dates <- sort(unique(c(solo.dates, double.dates)))
    
    # Apply model_fun to each ply
    runs <- bind_rows(
      get_runs(data.plys$odd.date.group),
      get_runs(data.plys$even.date.group)) %>%
      filter(date %in% unique.dates) %>%
      arrange(date) %>%
      mutate(date = as.Date(date, tz = tz(data$solar.time)))
    hour <- odd.date.group <- even.date.group <- ".dplyr.var"
    data_plys <- data.plys %>%
      select(-date, -hour, -odd.date.group, -even.date.group)
    out_list <- lapply(seq_along(runs$date), function(dt){
      run <- runs[dt,]
      data_ply <- data_plys[run$start:run$end,]
      data_daily_ply <-
        if(!is.null(data_daily)){
          data_daily[match(run$date, data_daily$date),]
        } else{
          NULL
        }
      ply_date <- run$date
      
      out <- model_fun(
        data_ply = data_ply, data_daily_ply = data_daily_ply,
        day_start = day_start, day_end = day_end, ply_date = ply_date,
        ply_validity = ply_validity, timestep_days = timestep_days,
        ...)
      # attach date column to output
      if(is.null(out) || nrow(out) == 0) NULL else data.frame(date=ply_date, out)
    })
  }
  
  # for when out_list comes back with nothing
  if(length(out_list) == 0 || all(sapply(out_list, is.null))){
    out_list <- tryCatch(
      list(
        model_fun(
          data_ply = data[c(),], data_daily_ply = data_daily[c(),],
          day_start = day_start, day_end = day_end, ply_date = as.Date(NA),
          ply_validity = NA, timestep_days = NA,
          ...) %>%
          mutate(date = as.Date(NA)) %>%
          select(date, everything())),
      error=function(e){
        data.frame(date = as.Date(NA))[c(),]
      }
    )
  }
  
  # Combine the 1-row dfs
  if(length(out_list) == 0){
    
  } else{
    example_choice <- 
      sapply(out_list, function(out){if(is.null(out)) 0 else ncol(out) }) %>%
      which.max()
    out_example <- out_list[[example_choice]][FALSE,]
    bind_rows(c(list(out_example), out_list)) %>% as.data.frame()
  }
  
}