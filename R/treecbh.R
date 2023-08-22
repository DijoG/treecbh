#' Function for stopping a program without error message.
stop_noerr <- function() {
  noerr = options(show.error.messages = FALSE)
  on.exit(options(noerr))
  stop()
}
#' Function for extraction point clouds to individual tree segments.
#' @param lasFILE las file of forest
#' @param multiPOLY sf multipolygon, individual tree segments
#' @param normalize logical, normalization based on CSF-classified ground points using the lidR::knnidw(), if FALSE "flat normalization" is performed
#' @param output_dir string, path to output directory
#' @param FEATURE character, attribute name to extract from multiPOLY to store in las
#' @return list of las files (point clouds of individual tree segments)
#' @export
get_3DTREE <- function(lasFILE, multiPOLY, normalize = TRUE, output_dir, FEATURE) {

  llas <- list()
  lground <- list()

  for (i in 1:(multiPOLY %>% nrow)) {
    llas[[i]] = clip_roi(lasFILE, multiPOLY[i,])
    if (normalize) {
      incsf = csf(sloop_smooth = TRUE, class_threshold = 1, cloth_resolution = 1, time_step = 1)
      llas[[i]] =
        llas[[i]] %>%
        filter_duplicates(.) %>%
        classify_ground(., incsf) %>%
        normalize_height(., knnidw()) %>%
        filter_poi(., Z > 0)
    } else {
      ZZ = llas[[i]]@data %>%
        filter(Classification == 2) %>%
        pull(Z)
      lground[[i]] = tibble(minZ = min(ZZ),
                            maxZ = max(ZZ),
                            meanZ = mean(ZZ))
      llas[[i]] = add_lasattribute(llas[[i]],
                                   llas[[i]]@data %>%
                                     mutate(Zn = Z - lground[[i]]$meanZ),
                                   name = "Zn",
                                   desc = "Zn") %>%
        filter_poi(., Zn >= 0) %>%
        filter_poi(., Classification != 2)
    }
    llas[[i]] = add_lasattribute(llas[[i]],
                                 multiPOLY[i, ] %>%
                                   pull(rlang::quo_squash(FEATURE)),
                                 name = rlang::quo_squash(FEATURE),
                                 desc = rlang::quo_squash(FEATURE))

    #> write out its las
    setwd(output_dir)
    writeLAS(llas[[i]], str_c("tree_0", i, ".las"))
  }
  return(llas)
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
#'
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
get_SEG <- function(list_LAS_char,
                    outdir1,
                    outdir2,
                    min_RANGE = 5,
                    min_POINT = .2,
                    K1 = 10, L1 = 1, DEC_R1 = .1,                   # First stage cut-pursuit parameters (treeiso)
                    K2 = 20, L2 = 20, MAX_GAP = .5, DEC_R2 = .1,    # Second stage cut-pursuit parameters (treeiso)
                    VER_O_W = .3, RHO = .5,                         # Final stage (treeiso)
                    cc_dir) {

  laso = list()
  for (tree in 1:length(list_LAS_char)) {
    message(crayon::silver(str_c("_______", tree, "________\n")))

    CC(cc_TREEiso(list_LAS_char[tree],
                  K1 = K1, L1 = L1, DEC_R1 = DEC_R1,
                  K2 = K2, L2 = L2, MAX_GAP = MAX_GAP, DEC_R2 = DEC_R2,
                  VER_O_W = VER_O_W, RHO = RHO,
                  output_dir = str_c(outdir1, basename(list_LAS_char[tree]))),
       cc_dir = cc_dir)

    las = readLAS(str_c(outdir1, basename(list_LAS_char[tree])))

    #> Remove points below 0.2 >
    lasl =
      las %>%
      filter_poi(., Z > min_POINT)

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
    if (length(msiz) <= 3) {

      #>
      if (length(msiz) == 0) {
        m =
          lasl@data %>%
          group_by(final_segs) %>%
          summarise(n = n(),
                    min_Z = min(Z),
                    max_Z = max(Z),
                    mean_Z = mean(Z),
                    range_Z = max(Z)-min(Z)) %>% # range_Z
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
                      range_Z = max(Z)-min(Z)) %>% # range_Z
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
                      range_Z = max(Z)-min(Z)) %>% # range_Z
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
      writeLAS(lasout, str_c(outdir2, basename(list_LAS_char[tree])))
      laso[[tree]] = lasout

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
        mm = mm
      } else {

        if (!"final_segs" %in% names (lasl@data)) {
          m =
            lasl@data %>%
            group_by(intermediate_segs) %>%
            summarise(n = n(),
                      min_Z = min(Z),
                      max_Z = max(Z),
                      mean_Z = mean(Z),
                      range_Z = max(Z)-min(Z)) %>% # range_Z
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
                      range_Z = max(Z)-min(Z)) %>% # range_Z
            data.frame()
          ms = sort(m$min_Z, index.return = TRUE)
          msi = lapply(ms, `[`, ms$x %in% head(unique(ms$x), nrow(m)))

          mm =
            m %>%
            filter(min_Z %in% msi$x)

          i_s = mm[mm$n == max(mm$n),]$final_segs

          lasout = lasl %>% filter_poi(., final_segs == i_s)

        }
        writeLAS(lasout, str_c(outdir2, basename(list_LAS_char[tree])))
        laso[[tree]] = lasout
      }

      #> Expectation rate (r) of BCH in relation to tree height, from bottom to top (default: .5)
      if (any(diff(mm$min_Z) > max(m$max_Z)/(exp(1-.5)*2.25))) { # or simply 2.25, mean(mm$range_Z)
        m =
          lasl@data %>%
          group_by(final_segs) %>%
          summarise(n = n(),
                    min_Z = min(Z),
                    max_Z = max(Z),
                    mean_Z = mean(Z),
                    range_Z = max(Z)-min(Z)) %>% # range_Z
          data.frame()
        ms = sort(m$min_Z, index.return = TRUE)
        msi = lapply(ms, `[`, ms$x %in% head(unique(ms$x), nrow(m)))

        mm =
          m %>%
          filter(range_Z ==  m$range_Z %>% max)

        i_s = mm$final_segs

        lasout = lasl %>% filter_poi(., final_segs == i_s)

      } else {

        #> Eliminate possible duplicates of min_Z
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

        #> Find the max range_Z, max_N and min_Z indices
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
          i_s = i_sR # or i_sM 3?
        }

        lasout = lasl %>% filter_poi(., intermediate_segs == i_s)
      }
      writeLAS(lasout, str_c(outdir2, basename(list_LAS_char[tree])))
      laso[[tree]] = lasout
    }
  }
  message(crayon::silver(str_c("_______ Done ________")))
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

#' Function for 2D kernel method including the decision scheme, used in get_CBH().
#' @param hist_DAT centered bins of counts on the output by	the vertical cross-sectional K-means clustering
#' @return numeric value of CBH
#' @export
get_CANOPYBH <- function(hist_DAT) {

  #hist_DAT =datt_histo
  #> Cluster kernel(s), finding how many kernels are >
  gp =
    ggplot(hist_DAT,
           aes(y = y, x = count)) +
    stat_density2d_filled(bins = 3) +
    geom_point()
  gpb = ggplot_build(gp)
  gpbdf =
    gpb[[1]][[1]]

  n_n =
    gpbdf %>%
    filter(piece == max(gpbdf$piece)) %>%
    pull(subgroup) %>%
    unique()

  if (gpbdf$x %>% max <= 10) {
    stop(return(height = gpbdf$y %>% max + 1))
  }

  #> One kernel >
  if (length(n_n) == 1) {
    z_r =
      gpbdf %>%
      filter(piece == max(gpbdf$piece)) %>%
      pull(y) %>%
      range() %>%
      diff()

    z_mi =
      gpbdf %>%
      filter(piece == max(gpbdf$piece)) %>%
      pull(y) %>%
      min()

    z_m =
      gpbdf %>%
      filter(piece == max(gpbdf$piece)) %>%
      pull(y) %>%
      max()

    if (z_r > gpbdf$y %>% range %>% diff *.5 &
        z_mi < gpbdf$y %>% range %>% max * .2 + (gpbdf$y %>% range %>% min)) {
      gpbdf_e =
        hist_DAT %>%
        filter(y >= z_mi & y <= z_m)
      max_n_ze = gpbdf_e[,"count"] %>% max()

      hist_ne = gpbdf_e[,"count"]

      rollone = gpbdf_e[,"y"] %>%
        round(., 3)

      heighte = rollone[which(hist_ne == max_n_ze)]
      if (length(heighte) > 1) {
        heighte = max(heighte)
      }

    } else {
      # > Cheking 2D space below and above kernel >
      if (gpbdf$y %>% max - z_m >= z_mi - gpbdf$y %>% min) {
        gpbdf_e =
          hist_DAT %>%
          filter(y >= z_m)
      } else {
        gpbdf_e =
          hist_DAT %>%
          filter(y <= z_mi)
      }

      max_n_ze = gpbdf_e[,"count"] %>% max()

      hist_ne = gpbdf_e[,"count"]

      rollone = gpbdf_e[,"y"] %>%
        round(., 3)

      heighte = rollone[which(hist_ne == max_n_ze)]
      if (length(heighte) > 1) {
        heighte = max(heighte)
      }
    }
  }

  max_n_z = hist_DAT[,"count"] %>% max()

  hist_n = hist_DAT[,"count"]

  rollon = hist_DAT[,"y"] %>%
    round(., 3)

  #> More than one kernel >
  if (length(n_n) > 1) {

    zrl = c()
    zmil = c()
    zml = c()
    kernel = c()
    for (i in 1:length(n_n)) {
      z_r_ =
        gpbdf %>%
        filter(piece == max(gpbdf$piece)) %>%
        filter(subgroup == i) %>%
        pull(y) %>%
        range() %>%
        diff()

      z_mi_ =
        gpbdf %>%
        filter(piece == max(gpbdf$piece)) %>%
        filter(subgroup == i) %>%
        pull(y) %>%
        min()

      z_m_ =
        gpbdf %>%
        filter(piece == max(gpbdf$piece)) %>%
        filter(subgroup == i) %>%
        pull(y) %>%
        max()
      zrl = c(zrl, z_r_)
      zmil = c(zmil, z_mi_)
      zml = c(zml, z_m_)
      kernel = c(kernel, i)
    }
    dfnz =
      tibble(z_r = zrl,
             z_min = zmil,
             z_max = zml,
             kernel = kernel)

    heights = c()
    for (i in 1:nrow(dfnz)) {
      gpbdf_e =
        hist_DAT %>%
        filter(y >= dfnz$z_min[i] & y <= dfnz$z_max[i])

      max_n_ze = gpbdf_e[,"count"] %>% max()

      hist_ne = gpbdf_e[,"count"]

      rollone = gpbdf_e[,"y"] %>%
        round(., 3)

      heighte = rollone[which(hist_ne == max_n_ze)]
      if (length(heighte) > 1) {
        heighte = heighte %>% min()
      }

      heights = c(heights, heighte)
    }
    dfnz =
      dfnz %>%
      mutate(heights = heights)
  }

  #> Decision >
  if (any(hist_n > 3)) {

    # > One kernel >
    if (length(n_n) == 1) {
      hist_n = hist_n[1:which(hist_n == max_n_z)]

      rollon = rollon[1:which(hist_n == max_n_z)]
      if (length(rollon) < 3) {
        height1 = z_m
      } else {
        # > Max count height along vertical profile >
        height1 = rollon[which(hist_n == max_n_z)]
      }

      if (length(height1) > 1) {
        height1 = max(height1)
      }

      if (height1 != heighte &
          z_r > gpbdf$y %>% range %>% diff *.5) {
        #> Inside/outside kernel max count height >
        height = heighte
      } else {
        #> Max count height of vertical profile >
        height = height1
      }
    } else {

      #> In case of more than one kernel >
      if (dfnz[dfnz$z_r == max(dfnz$z_r),]$z_min < gpbdf$y %>% range %>% diff *.4 &
          dfnz[dfnz$z_r == max(dfnz$z_r),]$z_r > gpbdf$y %>% range %>% diff *.5) {
        #> Max count height of longer kernel >
        height = dfnz[dfnz$z_r == max(dfnz$z_r),]$heights
      } else {
        #> Max count height of smaller kernel >
        height = dfnz$heights %>% min()
      }
    }
  } else {
    height = rollon[length(rollon)]
  }

  if (length(height) > 1) {
    height = height %>% min
  }

  return(height-(2*.2))
}

#' MAIN FUNCTION of treecbh, detecting CBH and deriving numerous metrics.
#' @param list_LAS_char character, list of las files
#' @param min_RANGE numeric, minimum height range (m, default = 5) of 3D tree segment employed during the process of within-segment tree isolation
#' @param min_POINT numeric, minimum height of points to eliminate forest floor and low vegetation (default = 0.2 m)
#' @param min_H_scale numeric, height scaler (m, default = .13), controlling understory removal
#' @param branch_WIDTH numeric, assumed CBH branch width (m, default = 0.2), controlling bin width for counting points
#' @param cross_WIDTH numeric, width of cross-section (m, default = 5)
#' @param cbh_ONLY numeric, options for executing: 1~treeiso and cbh, 2~only treeiso, 3~only cbh detection (default = 1, meaning tree isolation and CBH detection are active)
#' @param kM logical, interactive K-means cluster k tuning, activated if 'kM' = TRUE (default = FALSE)
#' @param method character, optional additional attribute (default = NULL)
#' @param outdir1 string, path to output directory of treeiso segment results
#' @param outdir2 string, path to output directory of filtered segments (intermediate_segs and final_segs)
#' @param K1,L1,DEC_R1 first stage cut-pursuit parameters (treeiso), default values as indicated
#' @param K2-L2,MAX_GAP,DEC_R2 second stage cut-pursuit parameters (treeiso), default values as indicated
#' @param VER_O_W,RHO final stage treeiso parameters, default values as indicated
#' @param cc_dir string, path to CloudCompare.exe
#' @return tibble (Z_max, Z_mean, Z_sd, Z_N_points, N_points, CBH, Hull_area, Del_vol, Cube_vol, Sphere_vol and treeID)
#' @export
get_CBH <- function(list_LAS_char,
                    min_RANGE = 5,
                    min_POINT = .2,
                    min_H_scale = .13,
                    branch_WIDTH = .2,
                    cross_WIDTH = 5,
                    cbh_ONLY = 1,
                    kM = FALSE,
                    # Interactive K-means cluster k tuning, activated if 'kM' = TRUE
                    # predicted k will be shown (plot), 'Do you accept k?'
                    # if you enter: "y", "Y", "yes", "YES", "ye", "YE" ~ accepted
                    # if you enter anything else:
                    # 'Enter k:' will ask you for a number to provide for k
                    method = NULL,
                    outdir1,
                    outdir2,
                    K1 = 10, L1 = 1, DEC_R1 = .1,
                    K2 = 20, L2 = 20, MAX_GAP = .5, DEC_R2 = .1,
                    VER_O_W = .3, RHO = .5,
                    cc_dir) {


  #> Possible errors >
  if (!min_H_scale %in% seq(.13,.25,.01)) {
    stop(crayon::magenta("Parameter 'min_H_scale' accepts values from .13 to .25"))
  }
  if (branch_WIDTH < 0) {
    stop(crayon::magenta("Parameter 'branch_WIDTH' must be positiv"))
  }
  if (cross_WIDTH < 4) {
    stop(crayon::magenta("Parameter 'cross_WIDTH' must be minimum 4"))
  }
  if (!cbh_ONLY %in% 1:3) {
    stop(crayon::magenta("Parameter 'cbh_ONLY' accepts 1, 2 and 3"))
  }
  if (K1 < 3 | K1 > 50) {
    stop(crayon::magenta("Parameter 'K1' accepts values btw. 3 and 50"))
  }
  if (K2 < 5 | K2 > 50) {
    stop(crayon::magenta("Parameter 'K2' accepts values btw. 5 and 50"))
  }
  if (L1 < .1 | L1 > 40) {
    stop(crayon::magenta("Regularizing parameter 'L1' accepts values btw. .1 and 40"))
  }
  if (L2 < 5 | L2 > 40) {
    stop(crayon::magenta("Regularizing parameter 'L2' accepts values btw. 5 and 40"))
  }
  if (DEC_R1 != .1) {
    stop(crayon::magenta("Parameter 'DEC_R1' must be .1"))
  }
  if (DEC_R2 != .1) {
    stop(crayon::magenta("Parameter 'DEC_R2' must be .1"))
  }
  if (MAX_GAP < .5 | MAX_GAP > 5) {
    stop(crayon::magenta("Parameter 'MAX_GAP accepts values btw. .5 and 5"))
  }
  if (VER_O_W < 0 | VER_O_W > 1.1) {
    stop(crayon::magenta("Parameter 'VER_O_W' accepts values btw. 0 and 1.1"))
  }
  if (RHO < 0 | RHO > 2) {
    stop(crayon::magenta("Parameter 'RHO' accepts values btw. 0 and 2"))
  }

  #> outdir strings
  if (str_sub(outdir1, nchar(outdir1), nchar(outdir1)) != "/") {
    outdir1 = str_c(outdir1, "/")
  }
  if (str_sub(outdir2, nchar(outdir2), nchar(outdir2)) != "/") {
    outdir2 = str_c(outdir2, "/")
  }

  #> ensure order
  list_LAS_char = gtools::mixedsort(list_LAS_char)

  if (cbh_ONLY %in% 1:2) {
    #> 3D tree decomposition (segmentation) >
    get_SEG(list_LAS_char,
            outdir1,
            outdir2,
            min_RANGE = min_RANGE,
            min_POINT = min_POINT,
            K1 = K1, L1 = L1, DEC_R1 = DEC_R1,
            K2 = K2, L2 = L2, MAX_GAP = MAX_GAP, DEC_R2 = DEC_R2,
            VER_O_W = VER_O_W, RHO = RHO,
            cc_dir = cc_dir)
    if (cbh_ONLY == 2) {
      message(crayon::green(str_c("_______ Done ________")))
      stop_noerr()
    }
  }

  if (cbh_ONLY %in% c(1,3)) {

    list_LASS = list.files(outdir2, pattern = ".las", full.names = T) %>%
      gtools::mixedsort()

    need = map_chr(list_LAS_char, ~str_split_1(.x, "/")[str_split_1(outdir2, "/") %>% length])
    list_LASS = list_LASS[grep(str_c(need, collapse = "|"), list_LASS)]

    metrics = list()
    for (tree in 1:length(list_LASS)) {
      message(crayon::green(str_c("_______", basename(list_LASS)[tree], "________")))

      #> Original las >
      laso = readLAS(list_LAS_char[tree])

      #> Segmented las >
      lass = readLAS(list_LASS[tree])

      #> min H ~ segmented las >
      m_H = (lass@data$Z %>% quantile(., .1) %>% as.vector) ^min_H_scale *2

      #> Horizontal cross section ~ segmented las >
      crosss = get_CROSS(lass, cross_WIDTH = cross_WIDTH)

      #> Denstiy (height ~ Z) on original las >
      dens =
        laso@data %>%
        ggplot(aes(y = Z)) +
        geom_histogram(binwidth = .2, center = 1)
      dat_dens =
        ggplot_build(dens)
      dat_dens =
        dat_dens[[1]][[1]]

      #> Removing ground points (again, just in case) and preparing segmented data for clustering >
      df = crosss@data %>%
        select(X, Z) %>%
        filter(Z > min_POINT)

      #> K-means clustering
      #> Prediction strength of a clustering:
      #> Tibshirani, R. and Walther, G. (2005) Cluster Validation by Prediction Strength, Journal of Computational and Graphical Statistics, 14, 511-528.
      set.seed(123)
      opk = fpc::prediction.strength(scale(df), 2, 10, 30)
      opkk = map_dbl(opk$predcorr, ~mean(.x)) %>% replace_na(., 0.1)
      k = ifelse(nrow(df) < 460, 1, ifelse(df$Z %>% min > 4, which.max(opkk) + 1, which.max(opkk)))
      k = ifelse(nrow(df) > 5000, which.max(opkk) + 2, k)
      k = ifelse(nrow(df) > 9000, 1, k)

      km = kmeans(scale(df), k, nstart = 25)

      if (kM == FALSE) {
        k = k
      } else {
        kmp = factoextra::fviz_cluster(km, data = df,
                                       palette = c("steelblue", "gold", "limegreen","grey","deeppink","forestgreen","grey45","steelblue","yellow"),
                                       geom = "point",
                                       ellipse.type = "convex",
                                       ggtheme = theme_bw()) +
          scale_y_continuous(breaks = NULL) +
          scale_x_continuous(breaks = NULL) +
          theme(axis.title.x = element_blank(),
                axis.text.x = element_blank(),
                plot.title = element_blank(),
                axis.ticks = element_blank(),
                axis.text.y = element_blank(),
                axis.title.y = element_blank(),
                legend.title = element_blank(),
                legend.position = c(.9,.25),
                legend.background = element_rect(fill = "transparent"))
        print(kmp)
        message(crayon::green(str_c("Suggested k is ", k, ".")))
        kM_ = readline(prompt = "Do you accept k? ")
        if (kM_ %in% c("y", "Y", "yes", "YES", "ye", "YE")) {
          km = km
        } else {
          kM__ = readline(prompt = "Enter k: ")
          kM__ = as.double(kM__)
          #> Update K-means >
          km = kmeans(scale(df), kM__, nstart = 25)
        }
      }

      dfk =
        df %>%
        mutate(km = km$cluster)

      dfkm =
        dfk %>% group_by(km) %>%
        summarise(mZ = mean(Z)) %>%
        bind_cols(km$centers)

      dfkk =
        dfk %>%
        filter(km == dfkm[which.min(dfkm$mZ),]$km)

      dfk_histo =
        dfkk %>%
        ggplot(aes(y = Z)) +
        geom_histogram(center = T, binwidth = branch_WIDTH) # Example
      datt_histo =
        ggplot_build(dfk_histo)
      datt_histo =
        datt_histo[[1]][[1]]

      datt_histo =
        datt_histo %>%
        filter(y > m_H)

      #> Finding CBH >
      cbh = get_CANOPYBH(datt_histo)

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

      #> Exclude ground and low vegetation from density (height ~ Z) of original las >
      dat_histos =
        dat_dens %>%
        filter(y >= m_H)

      #> Collect metrics >
      metrics[[tree]] =
        tibble(
          Z_max      = laso@data$Z %>% max,
          Z_mean     = laso@data$Z %>% mean,
          Z_sd       = laso@data$Z %>% sd,
          Z_N_points = dat_histos[dat_histos$count %>% which.max(),]$y,
          N_points   = dat_histos[dat_histos$count %>% which.max(),]$count,
          CBH        = cbh,
          Hull_area  = convex_hull$area,
          Del_vol    = delaunaj$areas %>% sum,
          Cube_vol   = round(nrow(vmv) * (0.2^3), 3),
          Sphere_vol = round(nrow(vmv) * (4/3)*(pi*(0.1^3)), 3),
          treeID     = tree)
    }

    metrics =
      metrics %>%
      bind_rows()

    if (!is.null(method)) {
      metrics =
        metrics %>%
        mutate(method = method)
    }
    message(crayon::green(str_c("_______ Done ________")))
    return(metrics)
  }
}

#' Function for 2D cross-sectional plot.
#' @param las las file
#' @param col character, color of points, default = "grey25"
#' @param cross_WIDTH numeric, width of cross-section (m, default = 5)
#' @return displays ggplot
#' @export
plot_CROSS <- function(las, col = "grey25", ylab, cross_WIDTH = 5) {
  p1 = c(min(las@data$X), mean(las@data$Y))
  p2 = c(max(las@data$X), mean(las@data$Y))
  data_clip = clip_transect(las, p1, p2, cross_WIDTH)
  PLOT = ggplot(data_clip@data, aes(X,Z)) +
    geom_point(size = .5, col = col) +
    coord_equal() +
    labs(y = ylab) +
    scale_y_continuous(breaks = NULL) +
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

