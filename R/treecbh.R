#' Function for stopping a program without error message.
stop_noerr <- function() {
  noerr = options(show.error.messages = FALSE)
  on.exit(options(noerr))
  stop()
}

#' Function for extraction point clouds to individual tree segments.
#' @param lasFILE las file of forest
#' @param multiPOLY sf multipolygon, individual tree segments
#' @param normalize logical, if TRUE normalization is performed on the ENTIRE dataset before segmentation
#' @param output_dir string, path to output directory
#' @param FEATURE character, attribute name and its values
#' @param RETURN logical, whether to return the list of las or not (default = FALSE)
#' @return list of las files (point clouds of individual tree segments)
#' @export
get_3DTREE <- function(lasFILE, multiPOLY, normalize = TRUE, output_dir, FEATURE = NULL, RETURN = FALSE) {

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Normalizing the dataset
  if (normalize) {
    message(crayon::blue("Normalizing entire point cloud..."))

    # Use conservative CSF parameters
    incsf = csf(sloop_smooth = FALSE, class_threshold = 0.5, cloth_resolution = 0.5, time_step = 0.65)

    tryCatch({
      lasFILE = lasFILE %>%
        filter_duplicates() %>%
        classify_ground(incsf) %>%
        normalize_height(knnidw()) %>%
        filter_poi(Z >= 0)

      message(crayon::green("Successfully normalized entire point cloud"))
    }, error = function(e) {
      message(crayon::red(paste("Normalization failed:", e$message)))
      stop("Cannot proceed without successful normalization")
    })
  }

  llas = list()

  for (i in 1:(nrow(multiPOLY))) {
    message(crayon::blue(paste("Processing tree", i, "of", nrow(multiPOLY))))

    # Clip the tree segment from the (already normalized) point cloud
    llas[[i]] = clip_roi(lasFILE, multiPOLY[i,])

    # Skip if no points in segment
    if (npoints(llas[[i]]) == 0) {
      message(crayon::yellow(paste("No points in tree segment", i, "- skipping")))
      next
    }

    # Add feature attribute if requested
    if (!is.null(FEATURE)) {
      feature_value = multiPOLY[i, ] %>% dplyr::pull(!!rlang::ensym(FEATURE))
      llas[[i]] = add_lasattribute(
        llas[[i]],
        x = feature_value,
        name = as.character(FEATURE),
        desc = as.character(FEATURE)
      )
    }

    # Write out the LAS file
    output_file = file.path(output_dir, paste0("tree_", sprintf("%03d", i), ".las"))
    writeLAS(llas[[i]], output_file)
  }

  if (RETURN) {
    return(llas)
  } else {
    message(crayon::green("_______ Done ________"))
  }
}

#' Function for executing CloudCompare plug-in, used in get_SEG(), taken from:
#' https://github.com/GreKro/cloudcompare.
#' @param cc_function function, cc_TREEiso()
#' @param cc_dir string, path to CloudCompare.exe
#' @return -
#' @export
CC <- function(cc_function, cc_dir) {
  CC_cmd = str_c(cc_dir, cc_function, sep = " ")
  for (f in 1:length(CC_cmd)) {
    system(command = CC_cmd[f])
  }
}

#' Function for grabbing treeiso plug-in employed in CloudCompare, used in CC(), taken from:
#' https://github.com/GreKro/cloudcompare and modified.
#' @param LAS_char string, path to las file including the name of las file
#' @param output_dir string, path to output directory
#' @param global_shift...no_timestamp CloudCompare parameters: https://www.cloudcompare.org/doc/wiki/index.php/Command_line_mode
#' @return las file
#' @export
cc_TREEiso <- function(LAS_char,
                       K1 = 5, L1 = 1.0, DEC_R1 = .05,                      # Cut-pursuit stage one (treeiso)
                       K2 = 20, L2 = 20, MAX_GAP = 2, DEC_R2 = 0.1,         # Cut-pursuit stage two (treeiso)
                       VER_O_W = .5, RHO = .5,                              # Final stage (treeiso)
                       output_dir,
                       global_shift = F,
                       global_shift_type = "AUTO",
                       filter_sf = F,
                       filter_value = c(0, 1),
                       c_export_fmt = "LAS",
                       c_ext = "las",
                       silent = T,
                       no_timestamp = T) {

  if(class(LAS_char) != "character" ) {
    stop(str_c("Wrong file is a ", class(LAS_char),", must be a character"))
  }

  if(class(filter_sf) != "logical" ) {
    stop("Wrong filter_sf must be a logical value - TRUE or FALSE")
  }

  if(length(filter_value) != 2) {
    stop("Wrong filter_value must be vector with 2 values for minimum and maximum treshold in filter_sf tool")
  }

  if(class(global_shift) != "logical" ) {
    stop("Wrong global_shift must be a logical value - TRUE or FALSE")
  }

  if(!(global_shift_type %in% c("AUTO", "FIRST")) && (class(global_shift_type) != "numeric" || length(global_shift_type) != 3)) {
    stop("Wrong global_shift_type must be character value : 'AUTO' or 'FIRST' (avaliable since CC v.2.11) or numeric vector with 3 values for dimensions x,y,z")
  }

  if(class(no_timestamp) != "logical" ) {
    stop("Wrong no_timestamp must be a logical value - TRUE or FALSE")
  }

  if(!c_export_fmt %in% c("LAS", "LAZ")) {
    stop(str_c('Wrong c_export_fmt: ',c_export_fmt,' (',class(c_export_fmt),'), must be character value, one of these : LAS, LAZ'))
  }

  if((!c_ext %in% c("LAS", "LAZ") && !c_ext %in% tolower(c("LAS","LAZ")))) {
    stop(str_c('Wrong c_ext: ', c_ext,' (',class(c_ext),'),
               must be character value, one of these :
               LAS, LAZ
               or
               las, laz'))

  }

  auto_save_off = "-AUTO_SAVE OFF"

  s1 = ifelse(silent, "-SILENT", "")

  s2 = ifelse(no_timestamp, "-NO_TIMESTAMP", "")

  s3 = str_c("-C_EXPORT_FMT", c_export_fmt, "-EXT", c_ext, sep = " ")

  s4 = ifelse(global_shift, str_c("-O -GLOBAL_SHIFT", global_shift_type, LAS_char, sep = " "), str_c("-O", LAS_char, sep = " "))

  s5 = str_c("-TREEISO", "-K1", K1, "-LAMBDA1", L1, "-DECIMATE_RESOLUTION1", DEC_R1,
             "-K2", K2, "-LAMBDA2", L2, "-MAX_GAP", MAX_GAP, "-DECIMATE_RESOLUTION2", DEC_R2,
             "-VERTICAL_OVERLAP_WEIGHT", VER_O_W, "-RHO", RHO, sep = " ")

  s6 = ifelse(filter_sf, str_c("-FILTER_SF", filter_value[1], filter_value[2], sep = " "), "")

  s7 = str_c("-SAVE_CLOUDS" , "FILE", output_dir, sep = " ")

  cc_function = str_c(s1, s2, auto_save_off, s3, s4, s5, s6, s7, sep = " ")

  return(cc_function)
}

#' Function for implementing within segment tree isolation using treeiso, used in get_CBH():
#' Xi, Z.; Hopkinson, C. 3D Graph-Based Individual-Tree Isolation (Treeiso) from Terrestrial Laser Scanning Point Clouds.
#' Remote Sens. 2022, 14, 6116. https://doi.org/10.3390/rs14236116
#' @param list_LAS_char character, list of las files
#' @param outdir1 string, path to output directory of tree isolation (treeiso) segments
#' @param outdir2 string, path to output directory of filtered segments (intermediate_segs and final_segs)
#' @param min_RANGE numeric, minimum height range (m) of 3D tree segment employed during the process of within-segment tree isolation, default = 5
#' @param min_POINT numeric, minimum height of points to eliminate forest floor and low vegetation (default = 0.2 m)
#' @param K1,L1,DEC_R1,K2,L2,MAX_GAP,DEC_R2,VER_O_W,RHO see function get_CBH()
#' @param cc_dir string, path to CloudCompare.exe
#' @return las file
#' @export
#' Function for implementing within segment tree isolation using treeiso, used in get_CBH():
#' Xi, Z.; Hopkinson, C. 3D Graph-Based Individual-Tree Isolation (Treeiso) from Terrestrial Laser Scanning Point Clouds.
#' Remote Sens. 2022, 14, 6116. https://doi.org/10.3390/rs14236116
#' @param list_LAS_char character, list of las files
#' @param outdir1 string, path to output directory of tree isolation (treeiso) segments
#' @param outdir2 string, path to output directory of filtered segments (intermediate_segs and final_segs)
#' @param min_RANGE numeric, minimum height range (m) of 3D tree segment employed during the process of within-segment tree isolation, default = 5
#' @param min_POINT numeric, minimum height of points to eliminate forest floor and low vegetation (default = 0.2 m)
#' @param K1,L1,DEC_R1,K2,L2,MAX_GAP,DEC_R2,VER_O_W,RHO see function get_CBH()
#' @param cc_dir string, path to CloudCompare.exe
#' @return las file
#' @export
get_SEG <- function(list_LAS_char,
                    outdir1,
                    outdir2,
                    min_RANGE = 5,
                    min_POINT = .2,
                    K1 = 10, L1 = 1, DEC_R1 = .1,                   # First stage cut-pursuit parameters (treeiso)
                    K2 = 20, L2 = 20, MAX_GAP = .5, DEC_R2 = .1,    # Second stage cut-pursuit parameters (treeiso)
                    VER_O_W = .3, RHO = .5,                         # Final stage (treeiso)
                    cc_dir) {

  #> Create output directories if they do not exist
  if (!dir.exists(outdir1)) {
    dir.create(outdir1, recursive = TRUE)
    message(crayon::green(paste("Created directory:", outdir1)))
  }
  if (!dir.exists(outdir2)) {
    dir.create(outdir2, recursive = TRUE)
    message(crayon::green(paste("Created directory:", outdir2)))
  }

  laso = list()
  processed_count = 0
  skipped_count = 0

  for (tree in 1:length(list_LAS_char)) {
    current_file <- list_LAS_char[tree]
    message(crayon::silver(str_c("_______", basename(current_file), "________\n")))

    # Define output file paths for cleanup
    output_file_path <- file.path(outdir1, basename(current_file))
    final_output_path <- file.path(outdir2, basename(current_file))

    # Track if current file was processed successfully
    file_processed <- FALSE

    #> CHECK 1: Skip if file does not exist
    if (!file.exists(current_file)) {
      message(crayon::yellow(paste("File does not exist:", current_file, "- skipping")))
      skipped_count <- skipped_count + 1
      next
    }

    #> CHECK 2: Check file size to avoid very small files
    file_size <- file.info(current_file)$size
    if (file_size < 1000) {       # Less than 1KB
      message(crayon::yellow(paste("File too small (", file_size, "bytes):", basename(current_file), "- skipping")))
      skipped_count <- skipped_count + 1
      next
    }

    #> CHECK 3: Read and check point count before processing
    las_check <- tryCatch({
      readLAS(current_file)
    }, error = function(e) {
      message(crayon::red(paste("Error reading file:", basename(current_file), "-", e$message)))
      return(NULL)
    })

    if (is.null(las_check)) {
      skipped_count <- skipped_count + 1
      next
    }

    #> CHECK 4: Skip if too few points
    point_count <- npoints(las_check)
    if (point_count < 20) {
      message(crayon::yellow(paste("Too few points (", point_count, ") in:", basename(current_file), "- skipping")))
      skipped_count <- skipped_count + 1
      next
    }

    #> If all checks pass, proceed with CloudCompare processing
    tryCatch({
      # Process with CloudCompare
      CC(cc_TREEiso(current_file,
                    K1 = K1, L1 = L1, DEC_R1 = DEC_R1,
                    K2 = K2, L2 = L2, MAX_GAP = MAX_GAP, DEC_R2 = DEC_R2,
                    VER_O_W = VER_O_W, RHO = RHO,
                    output_dir = output_file_path),
         cc_dir = cc_dir)

      # CHECK 5: Verify the output file was created
      if (!file.exists(output_file_path)) {
        stop("CloudCompare did not create output file")
      }

      # CHECK 6: Verify output file has points
      las_output <- tryCatch({
        readLAS(output_file_path)
      }, error = function(e) {
        message(crayon::red(paste("Error reading CloudCompare output:", e$message)))
        return(NULL)
      })

      if (is.null(las_output)) {
        # Clean up the invalid output file
        if (file.exists(output_file_path)) file.remove(output_file_path)
        stop("Could not read CloudCompare output")
      }

      if (npoints(las_output) == 0) {
        # Clean up the empty output file
        if (file.exists(output_file_path)) file.remove(output_file_path)
        stop("No points after CloudCompare processing")
      }

      #> Continue with the original processing logic
      las = las_output

      #> Remove points below min_POINT >
      lasl =
        las %>%
        filter_poi(., Z > min_POINT)

      #> Skip if no points left after filtering
      if (npoints(lasl) == 0) {
        # Clean up the output file since it became empty after filtering
        if (file.exists(output_file_path)) file.remove(output_file_path)
        stop(paste("No points above", min_POINT, "m after filtering"))
      }

      #> Calculate min_Z and range_Z >
      m =
        lasl@data %>%
        group_by(intermediate_segs) %>%
        summarise(n = n(),
                  min_Z = min(Z),
                  max_Z = max(Z),
                  mean_Z = mean(Z),
                  range_Z = max(Z)-min(Z)) %>%
        data.frame()

      #> Get Z and its indices >
      nr = ifelse(nrow(m) > 3, nrow(m), ifelse(nrow(m) == 1, 1, nrow(m)))
      ms = sort(m$min_Z, index.return = TRUE)
      msi = lapply(ms, `[`, ms$x %in% head(unique(ms$x), nr))

      #> Get Z >
      msiz = msi$x
      lasout <- NULL

      if (length(msiz) <= 3) {

        if (length(msiz) == 0) {
          m =
            lasl@data %>%
            group_by(final_segs) %>%
            summarise(n = n(),
                      min_Z = min(Z),
                      max_Z = max(Z),
                      mean_Z = mean(Z),
                      range_Z = max(Z)-min(Z)) %>%
            data.frame()
          nr = ifelse(nrow(m) > 3, 3, ifelse(nrow(m) == 1, 1, nrow(m)-1))
          ms = sort(m$min_Z, index.return = TRUE)
          msi = lapply(ms, `[`, ms$x %in% head(unique(ms$x), nr))

          mm =
            m %>%
            filter(range_Z ==  m$range_Z %>% max)

          i_s = mm$final_segs

          lasout = lasl %>% filter_poi(., final_segs == i_s)

        } else {
          if (!"final_segs" %in% names (lasl@data)) {
            m =
              lasl@data %>%
              group_by(intermediate_segs) %>%
              summarise(n = n(),
                        min_Z = min(Z),
                        max_Z = max(Z),
                        mean_Z = mean(Z),
                        range_Z = max(Z)-min(Z)) %>%
              data.frame()
            ms = sort(m$min_Z, index.return = TRUE)
            msi = lapply(ms, `[`, ms$x %in% head(unique(ms$x), nrow(m)))

            mm =
              m %>%
              filter(min_Z %in% msi$x)

            i_s = m[m$min_Z == min(mm$min_Z),]$intermediate_segs

            lasout = lasl %>% filter_poi(., intermediate_segs == i_s)

          } else {
            m =
              lasl@data %>%
              group_by(final_segs) %>%
              summarise(n = n(),
                        min_Z = min(Z),
                        max_Z = max(Z),
                        mean_Z = mean(Z),
                        range_Z = max(Z)-min(Z)) %>%
              data.frame()
            ms = sort(m$min_Z, index.return = TRUE)
            msi = lapply(ms, `[`, ms$x %in% head(unique(ms$x), nrow(m)))

            mm =
              m %>%
              filter(min_Z %in% msi$x)

            i_s = mm[mm$n== max(mm$n),]$final_segs

            lasout = lasl %>% filter_poi(., final_segs == i_s)

          }
        }

      } else {

        #> First 3 to 5 min_Z's!
        nr = ifelse(nrow(m) > 4, ifelse(nrow(m) > 5, 5, 4), 3)
        ms = sort(m$min_Z, index.return = TRUE)
        msi = lapply(ms, `[`, ms$x %in% head(unique(ms$x), nr))

        mm =
          m %>%
          filter(min_Z %in% msi$x) %>%
          arrange(min_Z)

        if (any(mm$range_Z > min_RANGE)) {
          # Continue with current mm
        } else {

          if (!"final_segs" %in% names (lasl@data)) {
            m =
              lasl@data %>%
              group_by(intermediate_segs) %>%
              summarise(n = n(),
                        min_Z = min(Z),
                        max_Z = max(Z),
                        mean_Z = mean(Z),
                        range_Z = max(Z)-min(Z)) %>%
              data.frame()
            ms = sort(m$min_Z, index.return = TRUE)
            msi = lapply(ms, `[`, ms$x %in% head(unique(ms$x), nrow(m)))

            mm =
              m %>%
              filter(min_Z %in% msi)

            i_s = mm[mm$min_Z == min(mm$min_Z),]$intermediate_segs

            lasout = lasl %>% filter_poi(., intermediate_segs == i_s)

          } else {
            m =
              lasl@data %>%
              group_by(final_segs) %>%
              summarise(n = n(),
                        min_Z = min(Z),
                        max_Z = max(Z),
                        mean_Z = mean(Z),
                        range_Z = max(Z)-min(Z)) %>%
              data.frame()
            ms = sort(m$min_Z, index.return = TRUE)
            msi = lapply(ms, `[`, ms$x %in% head(unique(ms$x), nrow(m)))

            mm =
              m %>%
              filter(min_Z %in% msi$x)

            i_s = mm[mm$n == max(mm$n),]$final_segs

            lasout = lasl %>% filter_poi(., final_segs == i_s)

          }
        }

        #> Expectation rate (r) of BCH in relation to tree height, from bottom to top (default: .5)
        if (any(diff(mm$min_Z) > max(m$max_Z)/(exp(1-.5)*2.25))) {
          m =
            lasl@data %>%
            group_by(final_segs) %>%
            summarise(n = n(),
                      min_Z = min(Z),
                      max_Z = max(Z),
                      mean_Z = mean(Z),
                      range_Z = max(Z)-min(Z)) %>%
            data.frame()
          ms = sort(m$min_Z, index.return = TRUE)
          msi = lapply(ms, `[`, ms$x %in% head(unique(ms$x), nrow(m)))

          mm =
            m %>%
            filter(range_Z ==  m$range_Z %>% max)

          i_s = mm$final_segs

          lasout = lasl %>% filter_poi(., final_segs == i_s)

        } else {

          # Eliminate possible duplicates of min_Z
          if (any(mm$min_Z %>% diff == 0)) {
            mm_z = mm$min_Z[duplicated(mm$min_Z)] %>% unique

            if (mm_z %>% length > 1) {

              mm_z_out = c()
              for (i in 1:length(mm_z)) {
                mm_z_o = mm %>%
                  filter(min_Z %in% mm_z[i])
                mm_z_ou = mm_z_o[mm_z_o$n < max(mm_z_o$n), ]$intermediate_segs
                mm_z_out = c(mm_z_out, mm_z_ou)
              }
            } else {
              mm_z_ou = mm[mm$min_Z == mm_z, ]
              mm_z_out = mm_z_ou[mm_z_ou$n < max(mm_z_ou$n),]$intermediate_segs
            }
            mm =
              mm %>%
              filter(!intermediate_segs %in% mm_z_out)
          }

          # Find the max range_Z, max_N and min_Z indices
          max_RANGE = mm$range_Z %>% max()
          max_N = mm$n %>% max()
          min_HEIGHT = mm$min_Z %>% min()

          if (mm[mm$min_Z == min_HEIGHT,]$range_Z < 2) {
            min_HEIGHT = mm$min_Z[2]
          }
          if ("final_segs" %in% names(mm)) {
            i_sN = mm[mm$n == max_N,]$final_segs
            i_sR = mm[mm$range_Z == max_RANGE,]$final_segs
            i_sM = mm[mm$min_Z == min_HEIGHT,]$final_segs
          } else {
            i_sN = mm[mm$n == max_N,]$intermediate_segs
            i_sR = mm[mm$range_Z == max_RANGE,]$intermediate_segs
            i_sM = mm[mm$min_Z == min_HEIGHT,]$intermediate_segs
          }
          if (length(i_sN) > 1) {
            if ("final_segs" %in% names(mm)) {
              i_sN = mm[mm$range_Z == m[m$n == max_N,]$range_Z %>% max(),]$final_segs
            } else {
              i_sN = mm[mm$range_Z == m[m$n == max_N,]$range_Z %>% max(),]$intermediate_segs
            }
          }

          if (i_sN == i_sR & i_sN == i_sM) {
            i_s = i_sR
          } else if (i_sN == i_sM & i_sM != i_sR) {
            i_s = i_sM
          } else if (i_sN == i_sR & i_sM != i_sR) {
            i_s = i_sN
          } else {
            i_s = i_sR
          }

          lasout = lasl %>% filter_poi(., intermediate_segs == i_s)
        }
      }

      #> FINAL CHECK: Only write and count if lasout has points
      if (!is.null(lasout) && npoints(lasout) > 0) {
        writeLAS(lasout, final_output_path)
        laso[[tree]] = lasout
        file_processed <- TRUE
        message(crayon::green(paste("Successfully processed:", basename(current_file))))
      } else {
        # Clean up both output files if final result is empty
        if (file.exists(output_file_path)) file.remove(output_file_path)
        if (file.exists(final_output_path)) file.remove(final_output_path)
        stop("No points in final output")
      }

    }, error = function(e) {
      message(crayon::red(paste("Processing failed for", basename(current_file), ":", e$message)))

      #> COMPREHENSIVE CLEANUP: Remove any output files that were created during failed processing
      if (file.exists(output_file_path)) {
        file.remove(output_file_path)
        message(crayon::yellow(paste("Cleaned up failed output:", basename(output_file_path))))
      }
      if (file.exists(final_output_path)) {
        file.remove(final_output_path)
        message(crayon::yellow(paste("Cleaned up failed output:", basename(final_output_path))))
      }
    })

    #> Update counts only once per file
    if (file_processed) {
      processed_count <- processed_count + 1
    } else {
      skipped_count <- skipped_count + 1
    }
  }

  #> Summary message - VERIFY COUNTS WITH ACTUAL FILES >
  actual_files_outdir1 <- length(list.files(outdir1, pattern = "\\.las$"))
  actual_files_outdir2 <- length(list.files(outdir2, pattern = "\\.las$"))

  message(crayon::green(str_c("Successfully processed: ", processed_count, " files")))
  message(crayon::yellow(str_c("Skipped: ", skipped_count, " files")))
  message(crayon::blue(str_c("Actual files in outdir1: ", actual_files_outdir1)))
  message(crayon::blue(str_c("Actual files in outdir2: ", actual_files_outdir2)))
  message(crayon::silver(str_c("_______ Done ________")))

  #> Return only non-null results >
  laso = Filter(Negate(is.null), laso)
  return(laso)
}

#' Function for getting 2D cross-section from 3D point cloud, used in get_CBH().
#' @param las las file
#' @param cross_WIDTH numeric, width of cross-section, metric, default = 5
#' @return las file
#' @export
get_CROSS <- function(las, cross_WIDTH = 5) {
  p1 = c(min(las@data$X), mean(las@data$Y))
  p2 = c(max(las@data$X), mean(las@data$Y))
  las_hcross = clip_transect(las, p1, p2, cross_WIDTH)
  return(las_hcross)
}

#' Function for percentile-based canopy base height detection, used in get_CBH().
#' @param hist_DAT centered bins of counts on the output by the vertical cross-sectional K-means clustering
#' @return numeric value of CBH
#' @export
get_CANOPYBH <- function(hist_DAT) {

  # Ensure histogram data validity
  if (is.null(hist_DAT) || nrow(hist_DAT) < 3 ||
      all(hist_DAT$count == 0) || all(is.na(hist_DAT$count))) {
    warning("Insufficient histogram data for CBH detection")
    return(NA)
  }

  # Ensure consitent column names
  if ("y" %in% names(hist_DAT)) {
    height_col <- "y"
  } else if ("x" %in% names(hist_DAT)) {
    height_col <- "x"
  } else if ("ymin" %in% names(hist_DAT)) {
    height_col <- "ymin"
  } else {
    warning("Cannot find height column in histogram data")
    return(NA)
  }

  # Filter out zero counts and ensure we have enough data
  hist_DAT = hist_DAT %>%
    filter(count > 0, !is.na(count), !is.na(!!sym(height_col)))

  if (nrow(hist_DAT) < 3) {
    warning("Not enough valid histogram bins for CBH detection")
    return(hist_DAT[[height_col]][which.max(hist_DAT$count)] - 0.4)
  }

  # SIMPLE FALLBACK METHOD - No kernel density, but MORE ACCURATE PERCENTILES!
  tryCatch({
    # Method 1: Use 5th percentile of point distribution
    cumulative = cumsum(hist_DAT$count) / sum(hist_DAT$count)
    cbh_index = which(cumulative >= 0.05)[1]

    if (!is.na(cbh_index)) {
      cbh = hist_DAT[[height_col]][cbh_index] - 0.4
      return(max(cbh, min(hist_DAT[[height_col]])))  # Ensure CBH is not below minimum height
    }

    # Method 2: If percentile fails, use point with maximum count
    max_count_idx = which.max(hist_DAT$count)
    if (length(max_count_idx) > 0) {
      return(hist_DAT[[height_col]][max_count_idx] - 0.4)
    }

    # Method 3: Ultimate fallback - median height
    return(median(hist_DAT[[height_col]]) - 0.4)

  }, error = function(e) {
    warning("All CBH detection methods failed: ", e$message)
    return(NA)
  })
}

#' MAIN FUNCTION of treecbh, detecting CBH and deriving numerous metrics.
#' @param list_LAS_char character, list of input tree point cloud las files (outputted by get_3DTREE())
#' @param min_RANGE numeric, minimum height range (m, default = 5) of 3D tree segment employed during the process of within-segment tree isolation
#' @param min_POINT numeric, minimum height of points to eliminate forest floor and low vegetation (default = 0.2 m)
#' @param min_H_scale numeric, height scaler (m, default = .13), controlling understory removal
#' @param branch_WIDTH numeric, assumed CBH branch width (m, default = 0.2), controlling bin width for counting points
#' @param cross_WIDTH numeric, width of cross-section (m, default = 5)
#' @param cbh_ONLY numeric, options for executing: 1~treeiso and cbh, 2~only treeiso, 3~only cbh detection (default = 3, meaning CBH detection is active)
#' @param kM logical, automatic percentile-based (default 'kM' = TRUE) or interactive CBH tuning
#' @param cbh_BUFF logical, if 'kM' = FALSE, the vertical buffer around the CBH (CBH-cbh_BUFF and CBH+cbh_BUFF) displayed with dotted red lines (default = 0.5)
#' @param method character, optional additional attribute (default = NULL in combination with default 'cbh_ONLY = 3')
#' @param outdir1 string, path to output directory of treeiso segment results (default = NULL in combination with default 'cbh_ONLY = 3')
#' @param outdir2 string, path to output directory of filtered segments (intermediate_segs and final_segs) (default = NULL)
#' @param K1,L1,DEC_R1 first stage cut-pursuit parameters (treeiso), default values as indicated
#' @param K2,L2,MAX_GAP,DEC_R2 second stage cut-pursuit parameters (treeiso), default values as indicated
#' @param VER_O_W,RHO final stage treeiso parameters, default values as indicated
#' @param cc_dir string, path to CloudCompare.exe (default = NULL in combination with default 'cbh_ONLY = 3')
#' @param VOL logical, if TRUE returns Hull_area, Del_vol, Cube_vol and Sphere_vol too (default = FALSE)
#' @return tibble (Z_max, Z_mean, Z_sd, Z_N_points, N_points, CBH and treeID) plus (Hull_area, Del_vol, Cube_vol, Sphere_vol)
#' @export
get_CBH <- function(list_LAS_char,
                    min_RANGE = 5,
                    min_POINT = 0.2,
                    min_H_scale = 0.13,
                    branch_WIDTH = 0.2,
                    cross_WIDTH = 5,
                    cbh_ONLY = 3,
                    kM = TRUE,
                    cbh_BUFF = 0.5,
                    method = NULL,
                    outdir1 = NULL,
                    outdir2 = NULL,
                    K1 = 10,
                    L1 = 1,
                    DEC_R1 = 0.1,
                    K2 = 20,
                    L2 = 20,
                    MAX_GAP = 0.5,
                    DEC_R2 = 0.1,
                    VER_O_W = 0.3,
                    RHO = 0.5,
                    cc_dir = NULL,
                    VOL = FALSE) {

  #> Possible errors >
  if (!min_H_scale %in% seq(0.13, 0.25, 0.01)) {
    stop(crayon::magenta("Parameter 'min_H_scale' accepts values from .13 to .25"))
  }
  if (branch_WIDTH < 0) {
    stop(crayon::magenta("Parameter 'branch_WIDTH' must be positive"))
  }
  if (cross_WIDTH < 4) {
    stop(crayon::magenta("Parameter 'cross_WIDTH' must be minimum 4"))
  }
  if (!cbh_ONLY %in% 1:3) {
    stop(crayon::magenta("Parameter 'cbh_ONLY' accepts 1, 2 and 3"))
  }
  if (cbh_BUFF < 0.1 | cbh_BUFF > 1) {
    stop(crayon::magenta("Parameter 'cbh_BUFF' ranges from 0.1 to 1"))
  }
  if (K1 < 3 | K1 > 50) {
    stop(crayon::magenta("Parameter 'K1' accepts values btw. 3 and 50"))
  }
  if (K2 < 5 | K2 > 50) {
    stop(crayon::magenta("Parameter 'K2' accepts values btw. 5 and 50"))
  }
  if (L1 < 0.1 | L1 > 40) {
    stop(crayon::magenta("Regularizing parameter 'L1' accepts values btw. .1 and 40"))
  }
  if (L2 < 5 | L2 > 40) {
    stop(crayon::magenta("Regularizing parameter 'L2' accepts values btw. 5 and 40"))
  }
  if (DEC_R1 != 0.1) {
    stop(crayon::magenta("Parameter 'DEC_R1' must be .1"))
  }
  if (DEC_R2 != 0.1) {
    stop(crayon::magenta("Parameter 'DEC_R2' must be .1"))
  }
  if (MAX_GAP < 0.5 | MAX_GAP > 5) {
    stop(crayon::magenta("Parameter 'MAX_GAP accepts values btw. .5 and 5"))
  }
  if (VER_O_W < 0 | VER_O_W > 1.1) {
    stop(crayon::magenta("Parameter 'VER_O_W' accepts values btw. 0 and 1.1"))
  }
  if (RHO < 0 | RHO > 2) {
    stop(crayon::magenta("Parameter 'RHO' accepts values btw. 0 and 2"))
  }

  #> Enhanced directory validation
  if (cbh_ONLY %in% c(1, 2)) {
    if (is.null(outdir1) || is.null(outdir2) || is.null(cc_dir)) {
      stop(crayon::magenta("'outdir1', 'outdir2' and 'cc_dir' are required when cbh_ONLY = 1 or 2 (tree isolation)"))
    }
    # Only process path strings if directories are provided
    if (!is.null(outdir1) && str_sub(outdir1, -1) != "/") outdir1 = str_c(outdir1, "/")
    if (!is.null(outdir2) && str_sub(outdir2, -1) != "/") outdir2 = str_c(outdir2, "/")
  }

  #> Ensure order
  list_LAS_char = gtools::mixedsort(list_LAS_char)

  if (cbh_ONLY %in% 1:2) {
    #> 3D tree decomposition (segmentation) >
    get_SEG(list_LAS_char, outdir1, outdir2, min_RANGE = min_RANGE,
            min_POINT = min_POINT, K1 = K1, L1 = L1, DEC_R1 = DEC_R1,
            K2 = K2, L2 = L2, MAX_GAP = MAX_GAP, DEC_R2 = DEC_R2,
            VER_O_W = VER_O_W, RHO = RHO, cc_dir = cc_dir)
    if (cbh_ONLY == 2) {
      message(crayon::green(str_c("_______ Done ________")))
      stop_noerr()
    }
  }

  if (cbh_ONLY %in% c(1, 3)) {
    #> File selection logic for interactive CBH detection
    if (cbh_ONLY == 1) {
      # After tree isolation, use files from outdir2
      list_LASS = list.files(outdir2, pattern = "\\.las$", full.names = TRUE) %>%
        gtools::mixedsort()
    } else if (cbh_ONLY == 3 && !is.null(outdir2)) {
      # CBH-only mode with specified output directory (pre-segmented trees)
      list_LASS = list.files(outdir2, pattern = "\\.las$", full.names = TRUE) %>%
        gtools::mixedsort()
    } else {
      # CBH-only mode without outdir2 - use original files for interactive detection
      list_LASS = list_LAS_char
    }

    # Validate whether there are files to process
    if (length(list_LASS) == 0) {
      stop(crayon::magenta("No LAS files found for CBH detection"))
    }

    metrics = list()
    for (tree in 1:length(list_LASS)) {
      message(crayon::green(str_c("_______", basename(list_LASS)[tree], "________")))

      #> Handle different file sources for interactive mode
      if (cbh_ONLY == 3 && is.null(outdir2)) {
        # Interactive CBH detection on original files
        laso = readLAS(list_LASS[tree])
        lass = laso  # Use the original file for segmentation analysis
      } else {
        # Normal flow: original file and segmented file from tree isolation
        laso = readLAS(list_LAS_char[tree])
        lass = readLAS(list_LASS[tree])
      }

      #> min H ~ segmented las >
      m_H = (lass@data$Z %>% quantile(., .1) %>% as.vector) ^ min_H_scale * 2

      #> Horizontal cross section ~ segmented las >
      crosss = get_CROSS(lass, cross_WIDTH = cross_WIDTH)

      #> Horizontal cross section ~ original las >
      cross = get_CROSS(laso, cross_WIDTH = cross_WIDTH)

      #> Density (height ~ Z) on original las >
      dens =
        laso@data %>%
        ggplot(aes(y = Z)) +
        geom_histogram(binwidth = .2)
      dat_dens =
        ggplot_build(dens)
      dat_dens =
        dat_dens[[1]][[1]]

      #> Removing ground points (again, just in case) and preparing segmented data for clustering >
      df = crosss@data %>%
        select(X, Z) %>%
        filter(Z > min_POINT)

      #> Removing points under min_POINT and preparing original data >
      dfo = cross@data %>%
        select(X, Z) %>%
        filter(Z > min_POINT)

      #> K-means clustering
      #> Prediction strength of a clustering:
      #> Tibshirani, R. and Walther, G. (2005) Cluster Validation by Prediction Strength, Journal of Computational and Graphical Statistics, 14, 511-528.
      set.seed(123)
      opk = fpc::prediction.strength(scale(df), 2, 10, 30)
      opkk = map_dbl(opk$predcorr, ~mean(.x)) %>% replace_na(., 0.1)

      k = which.max(opkk)
      km = kmeans(scale(df), k, nstart = 25)

      dfk = df %>% mutate(km = km$cluster)
      dfkm = dfk %>% group_by(km) %>% summarise(mZ = mean(Z)) %>%
        bind_cols(km$centers)
      dfkk = dfk %>% filter(km == dfkm[which.min(dfkm$mZ),]$km)

      #> FIXED: Remove 'center' parameter and handle column names properly >
      dfk_histo = dfkk %>% ggplot(aes(y = Z)) + geom_histogram(binwidth = branch_WIDTH)

      datt_histo = ggplot_build(dfk_histo)
      datt_histo = datt_histo[[1]][[1]]

      #> FIXED: Handle different column names from ggplot_build >>
      if ("y" %in% names(datt_histo)) {
        datt_histo = datt_histo %>% filter(y > m_H)
      } else if ("x" %in% names(datt_histo)) {
        datt_histo = datt_histo %>% filter(x > m_H)
      } else if ("ymin" %in% names(datt_histo)) {
        datt_histo = datt_histo %>% filter(ymin > m_H)
      } else {
        # Use first numeric column as fallback
        numeric_cols = sapply(datt_histo, is.numeric)
        if (any(numeric_cols)) {
          height_col = names(datt_histo)[which(numeric_cols)[1]]
          datt_histo = datt_histo %>% filter(!!sym(height_col) > m_H)
        }
      }

      #> Finding CBH >
      cbh = get_CANOPYBH(datt_histo)

      #> Activation of CBH tuning >
      if (kM) {
        cbh = cbh
      } else {
        p = plot_CROSS(laso, ylab = "Height [m]") +
          geom_hline(yintercept = cbh, col = "firebrick3") +
          geom_hline(yintercept = cbh - cbh_BUFF, col = "firebrick2", linetype = "dotted") +
          geom_hline(yintercept = cbh + cbh_BUFF, col = "firebrick2", linetype = "dotted")
        print(p)

        message(crayon::green(str_c("Suggested CBH is ", round(cbh, 2), ".")))

        # FIXED: Simpler input handling without complex validation loops
        kM_ <- readline(prompt = "Do you accept this CBH? (y = yes, n = no): ")

        if (tolower(trimws(kM_)) %in% c("y", "yes", "ye", "")) {
          # User accepts the suggested CBH
          cbh = cbh
          message(crayon::green("Using suggested CBH."))
        } else {
          # User wants to enter custom CBH
          message("Enter your assumed CBH.")
          kM__ <- readline(prompt = "Assumed CBH: ")

          # Convert to numeric and validate
          kM__ <- as.numeric(kM__)

          if (is.na(kM__) || !is.numeric(kM__)) {
            message(crayon::red("Invalid input. Using suggested CBH."))
            # Keep the original cbh value
          } else if (kM__ < 0 || kM__ > max(laso@data$Z, na.rm = TRUE)) {
            message(crayon::red("CBH out of reasonable range. Using suggested CBH."))
            # Keep the original cbh value
          } else {
            # Valid custom CBH entered
            dfk_histog = dfo %>%
              filter(Z > (kM__ - cbh_BUFF) & Z < (kM__ + cbh_BUFF)) %>%
              ggplot(aes(y = Z)) +
              geom_histogram(binwidth = branch_WIDTH)

            datt_histog = ggplot_build(dfk_histog)
            datt_histog = datt_histog[[1]][[1]]

            #> Finding CBH based on assumed (tuned) CBH >
            cbh = get_CANOPYBH(datt_histog)
            message(crayon::green(paste("Using user-specified CBH:", round(cbh, 2))))
          }
        }
      }

      #> Exclude ground and low vegetation from density (height ~ Z) of original las >
      if ("y" %in% names(dat_dens)) {
        dat_histos = dat_dens %>% filter(y >= m_H)
      } else if ("x" %in% names(dat_dens)) {
        dat_histos = dat_dens %>% filter(x >= m_H)
      } else if ("ymin" %in% names(dat_dens)) {
        dat_histos = dat_dens %>% filter(ymin >= m_H)
      } else {
        # Use first numeric column as fallback
        numeric_cols = sapply(dat_dens, is.numeric)
        if (any(numeric_cols)) {
          height_col = names(dat_dens)[which(numeric_cols)[1]]
          dat_histos = dat_dens %>% filter(!!sym(height_col) >= m_H)
        } else {
          dat_histos = dat_dens
        }
      }

      #> Collect metrics >
      if (VOL) {
        #> Canopy ~ original las >
        newdata =
          laso@data %>%
          data.frame() %>%
          select(X, Y, Z) %>%
          filter(Z > cbh) %>%
          as.matrix()

        las_new =
          laso %>%
          filter_poi(Z > cbh)

        #> Delaunay hull area >
        convex_hull =
          geometry::convhulln(newdata, "FA")

        #> Delaunay hull volume >
        delaunaj =
          geometry::delaunayn(newdata, "Fa")

        #> Voxelizing canopy (0.2m) >
        vmv = get_VOXEL(las_new@data, .2)

        #> FIXED: Handle column names in dat_histos for metrics >
        if ("y" %in% names(dat_histos)) {
          max_count_height = dat_histos$y[which.max(dat_histos$count)]
          max_count_value = dat_histos$count[which.max(dat_histos$count)]
        } else if ("x" %in% names(dat_histos)) {
          max_count_height = dat_histos$x[which.max(dat_histos$count)]
          max_count_value = dat_histos$count[which.max(dat_histos$count)]
        } else {
          max_count_height = NA
          max_count_value = NA
        }

        metrics[[tree]] =
          tibble(
            Z_max      = laso@data$Z %>% max,
            Z_mean     = laso@data$Z %>% mean,
            Z_sd       = laso@data$Z %>% sd,
            Z_N_points = max_count_height,
            N_points   = max_count_value,
            CBH        = cbh,
            Hull_area  = convex_hull$area,
            Del_vol    = delaunaj$areas %>% sum,
            Cube_vol   = round(nrow(vmv) * (0.2^3), 3),
            Sphere_vol = round(nrow(vmv) * (4/3)*(pi*(0.1^3)), 3),
            treeID     = str_remove(basename(list_LASS)[tree], ".las") %>% str_remove(., "tree_"))
      } else {
        #> FIXED: Handle column names in dat_histos for metrics >
        if ("y" %in% names(dat_histos)) {
          max_count_height = dat_histos$y[which.max(dat_histos$count)]
          max_count_value = dat_histos$count[which.max(dat_histos$count)]
        } else if ("x" %in% names(dat_histos)) {
          max_count_height = dat_histos$x[which.max(dat_histos$count)]
          max_count_value = dat_histos$count[which.max(dat_histos$count)]
        } else {
          max_count_height = NA
          max_count_value = NA
        }

        metrics[[tree]] =
          tibble(
            Z_max      = laso@data$Z %>% max,
            Z_mean     = laso@data$Z %>% mean,
            Z_sd       = laso@data$Z %>% sd,
            Z_N_points = max_count_height,
            N_points   = max_count_value,
            CBH        = cbh,
            treeID     = str_remove(basename(list_LASS)[tree], ".las") %>% str_remove(., "tree_"))
      }
    }

    metrics = metrics %>% bind_rows()

    if (!is.null(method)) {
      metrics = metrics %>% mutate(method = method)
    }
    message(crayon::green(str_c("_______ Done ________")))
    return(metrics)
  }
}

#' Automated parameter configuration for get_CBH() based on point cloud density
#' @param input_dir string, path to directory containing LAS files
#' @param density_threshold numeric, point density threshold (points/m²) for automatic vs interactive mode (default = 20)
#' @return data.frame with recommended parameters for get_CBH
#' @export
get_PARAMS <- function(input_dir, density_threshold = 20) {

  # Input validation
  if (!dir.exists(input_dir)) {
    stop("Input directory does not exist: ", input_dir)
  }

  las_files = list.files(input_dir, pattern = "\\.las$", full.names = TRUE) %>%
    gtools::mixedsort()

  if (length(las_files) == 0) {
    stop("No LAS files found in: ", input_dir)
  }

  # Calculate densities with error handling
  densities = sapply(las_files, function(file) {
    tryCatch({
      las = readLAS(file)
      lidR::density(las)
    }, error = function(e) {
      warning("Could not read file: ", basename(file), " - ", e$message)
      NA
    })
  })

  # Remove NA values and calculate statistics
  valid_densities = densities[!is.na(densities)]

  if (length(valid_densities) == 0) {
    stop("No valid LAS files could be processed")
  }

  mean_density = mean(valid_densities, na.rm = TRUE)
  min_density = min(valid_densities, na.rm = TRUE)
  max_density = max(valid_densities, na.rm = TRUE)

  # Determine optimal parameters
  if (mean_density < density_threshold) {
    # Low density: Use interactive mode on original files
    cbh_ONLY = 3
    kM = FALSE
    recommendation = "Low point density detected. Using interactive CBH detection on original files for manual verification."
  } else {
    # High density: Use automatic tree isolation + CBH
    cbh_ONLY = 1
    kM = TRUE
    recommendation = "High point density detected. Using automatic tree isolation and CBH detection."
  }

  # Create comprehensive results
  result = list(
    parameters = data.frame(
      cbh_ONLY = cbh_ONLY,
      kM = kM,
      density_threshold_used = density_threshold
    ),
    density_stats = data.frame(
      mean_density = round(mean_density, 2),
      min_density = round(min_density, 2),
      max_density = round(max_density, 2),
      n_files_processed = length(valid_densities),
      n_files_total = length(las_files)
    ),
    recommendation = recommendation
  )

  # Print summary for user
  message("=== get_PARAMS Analysis ===")
  message("Density statistics (points/m²):")
  message("  Mean: ", round(mean_density, 2))
  message("  Range: ", round(min_density, 2), " - ", round(max_density, 2))
  message("Recommendation: ", recommendation)
  message("Suggested parameters:")
  message("  cbh_ONLY = ", cbh_ONLY)
  message("  kM = ", kM)

  return(result)
}

#' Function for 2D cross-sectional plot used only in get_CBH().
#' @param las las file
#' @param col character, color of points, default = "grey25"
#' @param cross_WIDTH numeric, width of cross-section (m, default = 5)
#' @return displays ggplot
#' @export
plot_CROSS <- function (las, col = "grey25", ylab, cross_WIDTH = 5) {
  p1 = c(min(las@data$X), mean(las@data$Y))
  p2 = c(max(las@data$X), mean(las@data$Y))
  data_clip = clip_transect(las, p1, p2, cross_WIDTH)
  PLOT = ggplot(data_clip@data, aes(X, Z)) +
    geom_point(size = 0.5, col = col) +
    coord_equal() +
    labs(y = ylab, x = NULL) +
    scale_y_continuous(breaks = ceiling(min(data_clip$Z)):floor(max(data_clip$Z))) +
    scale_x_continuous(breaks = NULL)
  return(PLOT)
}

#' Helper functions from VOXR package:
#' Bastien Lecigne and others, Exploring trees in three dimensions:
#' VoxR, a novel voxel-based R package dedicated to analysing the complex arrangement of tree crowns,
#' Annals of Botany, Volume 121, Issue 4, 14 March 2018, Pages 589–601,
#' https://doi.org/10.1093/aob/mcx095
#'
#' Function for voxelizing point cloud.
#' @importFrom data.table ":="
#' @param data data.table, tibble, data.frame  as las@data
#' @param res_VOXEL numeric, voxel resolution
#' @param full.grid logical, if TRUE empty voxels in the tree bounding box are returned
#' @param message logical, if TRUE error and interactive messages are enabled
#' @return data table
#' @export
get_VOXEL <- function (data, res_VOXEL, full.grid, message) {
  X = Y = Z = npts = .N = . = `:=` = NULL
  if (!(data.table::is.data.table(data))) {
    data = data.table::data.table(data)
  }
  check = check_DAT(data, message = message)
  if (missing(res_VOXEL)) {
    stop("No voxel resolution 'res_VOXEL' provided")
  }
  else {
    if (!is.vector(res_VOXEL))
      stop("res_VOXEL must be a vector of length 1")
    if (!is.numeric(res_VOXEL))
      stop("res_VOXEL must be numeric")
    if (res_VOXEL <= 0)
      stop("res_VOXEL must be positive")
    if (length(res_VOXEL) > 1) {
      res_VOXEL = res_VOXEL[1]
      warning("res_VOXEL contains more than 1 element. Only the first was used")
    }
  }
  if (missing(full.grid))
    full.grid = FALSE
  data = check$data
  data[, `:=`(X = Rfast::Round(X/res_VOXEL) * res_VOXEL,
              Y = Rfast::Round(Y/res_VOXEL) * res_VOXEL,
              Z = Rfast::Round(Z/res_VOXEL) * res_VOXEL)]
  data = unique(data[, `:=`(npts, .N), by = .(X, Y, Z)])
  if (full.grid) {
    x_seq = seq(min(data$X), max(data$X), res_VOXEL)
    y_seq = seq(min(data$Y), max(data$Y), res_VOXEL)
    z_seq = seq(min(data$Z), max(data$Z), res_VOXEL)
    est_weight = round(length(x_seq) * length(y_seq) * length(z_seq) *
                         4 * 8/2^{20}/1024, 1)
    if (est_weight > 2) {
      cat(paste("Final data is estimated to be ", est_weight,
                "GB in size. Continue execution ? y or n"))
      test = readline()
      if (test != "y") {
        stop("Execution stoped.")
      }
    }
    empty = data.table::data.table(expand.grid(x_seq,
                                               y_seq,
                                               z_seq))
    data.table::setnames(empty, c("X", "Y", "Z"))
    empty[, `:=`(npts, 0)]
    data = dplyr::bind_rows(data, empty)
    data = data[, c(`:=`(npts, sum(npts)),
                    `:=`(Intensity, mean(Intensity))), keyby = .(X, Y, Z)]
  }
  if (check$dfr)
    data = as.data.frame(data)
  return(data)
}
#'
#' Function for checking input data used in get_VOXEL().
#' @param data data.table, tibble, data.frame  as las@data
#' @param message logical, if TRUE error and interactive messages are enabled
#' @return data table
#' @export
check_DAT <- function (data, message) {
  if (missing(message))
    message = FALSE
  if (!(is.data.frame(data)))
    stop("data must be a data.frame or a data.table")
  if (ncol(data) > 3) {
    if (message)
      print("NOTE: data contain more than 3 columns, three first and 'Intensity' are used")
  }
  dfr = FALSE
  if (!data.table::is.data.table(data))
    dfr = TRUE
  data = data.table::data.table(data[, c(1:3,5)])
  data.table::setnames(data, c("X", "Y", "Z", "Intensity"))
  data$Intensity = data$Intensity %>% as.numeric()
  if (!(all(sapply(data, class) == "numeric"))) {
    stop("All the fields of the data must be numeric")
  }
  if (any(is.na(data)) & message)
    warning("data contains missing values.")
  invisible(list(data = data, dfr = dfr))
}
