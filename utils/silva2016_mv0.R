# This script provides code to a modified lidR::silva2016() segmentation algorithm; 
# The function corrects 'max_cr_factor' input that reduces the max diameter for small plants, 
# it takes in 'hmax' and modifies it so that a correction factor can be applied to the 'max_cr_factor*hmax' product.

silva2016_mv0 <- function (chm, treetops, max_cr_factor = 0.6, exclusion = 0.3, 
          ID = "treeID", fun = NULL) 
{
  assert_is_a_number(max_cr_factor)
  assert_is_a_number(exclusion)
  assert_all_are_positive(max_cr_factor)
  assert_all_are_in_open_range(exclusion, 0, 1)
  treetops <- check_tree_tops(treetops, ID)
  chm <- lazyeval::uq(chm)
  treetops <- lazyeval::uq(treetops)
  max_cr_factor <- lazyeval::uq(max_cr_factor)
  exclusion <- lazyeval::uq(exclusion)
  ID <- lazyeval::uq(ID)
  f = function(bbox) {
    assert_is_valid_context(LIDRCONTEXTITS, "silva2016", 
                            null_allowed = TRUE)
    if (nrow(treetops) == 0L) {
      crown <- chm
      crown <- raster_set_values(crown, NA_integer_)
      return(crown)
    }
    if (raster_is_proxy(chm) & missing(bbox)) 
      stop("Cannot segment the trees from a raster stored on disk. Use segment_trees() or load the raster in memory", 
           call. = FALSE)
    res <- crop_special_its(treetops, chm, bbox)
    treetops <- res$treetops
    chm <- res$chm
    st_chm <- chm
    if (!is(chm, "stars")) 
      st_chm <- stars::st_as_stars(chm)
    if (nrow(treetops) == 0L) {
      warning("No tree can be used as seed", call. = FALSE)
      crown <- chm
      crown <- raster_set_values(crown, NA_integer_)
      return(crown)
    }
    chmdt <- raster_as_dataframe(chm, xy = FALSE, na.rm = TRUE)
    data.table::setDT(chmdt)
    ids <- treetops[[ID]]
    coords <- sf::st_coordinates(treetops)
    u <- C_knn(coords[, 1], coords[, 2], chmdt$X, chmdt$Y, 
               1L, getThread())
    id <- d <- hmax <- Z <- . <- X <- Y <- NULL
    chmdt[, `:=`(id, u$nn.idx[, 1])]
    chmdt[, `:=`(id, ids[id])]
    chmdt[, `:=`(d, u$nn.dist[, 1])]
    chmdt[, `:=`(hmax, max(Z)), by = id]

    chmdt <- chmdt[Z >= exclusion * hmax & d <= max_cr_factor * 
                     hmax * ifelse(is.function(fun), fun(hmax), 1), .(X, Y, id)]
    crown <- chm
    crown <- raster_set_values(crown, NA_integer_)
    cells <- raster_cell_from_xy(crown, chmdt$X, chmdt$Y)
    crown <- raster_set_values(crown, chmdt[["id"]], cells)
    return(crown)
  }
  f <- plugin_its(f, omp = TRUE, raster_based = TRUE)
  return(f)
}
# this line ensures that internally used packages in lidR are available for the modified function
environment(silva2016_mv0) <- asNamespace('lidR')
