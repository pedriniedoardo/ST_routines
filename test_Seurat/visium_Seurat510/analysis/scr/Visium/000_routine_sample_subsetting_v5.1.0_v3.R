# load sample dataset -----------------------------------------------------
# load a sample assay v5
test <- readRDS("../../data/misc/stxBrain.SeuratData_0.1.2_anterior1.rds")

# confirm the class of the object
class(test@assays$Spatial)

# I can pull the coordinates
GetTissueCoordinates(test)

Images(test)
test@images

# this fails: Spatial coordinates are no longer fetchable with FetchData
subset(test, slice1_x > 400 | slice1_y < 150, invert = TRUE)


# workaround --------------------------------------------------------------
# add the coordinates to the metadata to make it possible the subsetting
df_coordinates <- data.frame(barcode = test@images$anterior1$centroids@cells,
                             imagerow = test@images$anterior1$centroids@coords[, 1],
                             imagecol = test@images$anterior1$centroids@coords[, 2])

df_meta_full <- test@meta.data %>%
  rownames_to_column("barcode") %>%
  left_join(df_coordinates, by = "barcode") %>%
  column_to_rownames("barcode")


df_meta_full %>%
  ggplot(aes(y=imagerow, x=imagecol))+ geom_point() + scale_y_reverse()

# swap the metadata
test@meta.data <- df_meta_full

# now the subsetting works
test_subset <- subset(test,imagerow < 7000 & imagerow > 2500)
test_subset <- subset(test_subset,imagecol > 3500 & imagecol < 10000)

SpatialDimPlot(test_subset, crop = TRUE, label = TRUE)
