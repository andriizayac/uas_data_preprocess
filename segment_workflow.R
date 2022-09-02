pkgs <- c("lidR", "raster", "gstat", "sf", "sp", "viridisLite", "stars", "dplyr")
sapply(pkgs, require, character.only = TRUE)

# === helper stuff
cols <- gray.colors(255)
colsid <- sample(viridis(50, alpha = 1, direction = 1, begin = 0))

# ===
# ---read in the las catalog
path <- "D:/UAV_analysis/SodaBurn2021/"
# --- select a subset
tiles <- readLAScatalog(paste0(path, "tiles/")) 
ts <- catalog_select(tiles)
# --- retile if needed
# opt_output_files(ts) <- "C:/tmp/dtm_workflow_example/crp_regen3_2021_{ID}"
# opt_chunk_buffer(ts) <- 0
# opt_chunk_size(ts) <- 8
# newts <- catalog_retile(ts)

# bring dtm 
dtm <- raster(list.files(path, pattern = "20210615_DTM.*\\.tif", full.names = TRUE))

# normalize las
las <- readLAS(ts, select = "xyz")
# lasn <- normalize_height(las, dtm)
# lasn <- filter_poi(lasn, Z > -.2)

# import and crop the CHM
clist <- list.files(path, pattern = "_CHM_wgs84utmz11n.tif", full.names = TRUE)
chm <- sapply(clist, raster)
chms <- sapply(chm, crop, y = las@bbox)

# create roi (buffer = -.9m) to avoid half-captured plants
roi <- chms[[1]] %>% 
  st_bbox() %>% 
  st_as_sfc() %>% 
  st_buffer(-.9)

# find trees (each year separately)
tl <- sapply(chms, function(x) {
  find_trees(x, lmf(ws = .6, hmin = .05, shape = "circular"))
  })

# use a single layer for tree tops with the specified base layer
j <- 3
ttops <- list(tl[[j]], tl[[j]], tl[[j]])

# segment trees in CHM
cro <- sapply(1:3, function(x) {
  silva2016(chms[[x]], ttops[[x]], max_cr_factor = 1.5, exclusion = 0.25)()
  })

# === visualize
yr <- 3
plot(chms[[yr]], col = cols)
plot(cro[[yr]], add = TRUE, col = colsid)
plot(ttops[[yr]], add = TRUE, pch = 19, cex = .5, col = "black")
plot(roi, add = TRUE)


# --- subset to roi
tout <- lapply(1:3, function(x) {
  ttops[[x]] %>% 
  st_as_sf() %>% 
  st_crop(st_bbox(roi)) 
  })

cout <- lapply(1:3, function(x) {
  cro[[x]] %>% 
  st_as_stars() %>% 
  st_as_sf(as_points = FALSE, merge = TRUE) %>% 
  mutate(var = lengths(st_intersects(., st_buffer(tout[[x]], .1)) )) %>% 
  filter(var == 1) %>% 
  select(-var) 
  })

# === export
yrs <- c(2019:2021)
for(i in 1:length(yrs)){
  tout[[i]] %>% 
    st_write(paste0(path, "segment_output/ttops_", yrs[[i]], ".shp"), delete_layer = FALSE)
  
  cout[[i]] %>% 
    st_write(paste0(path, "segment_output/crowns_", yrs[[i]],".shp"), delete_layer = FALSE)
}


