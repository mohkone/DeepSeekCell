test_that("LLM JSON responses are parsed with trailing text", {
  response <- paste0(
    "```json\n",
    '{"annotations":[{"cluster":"Cluster1","cell_type":"T cell","confidence":95,',
    '"is_mixed":"false","tissue_consistency":"expected","reasoning":"CD3 markers"}]}',
    "\n```\nextra text"
  )

  parsed <- parse_annotation_response(response)

  expect_s3_class(parsed, "data.frame")
  expect_equal(nrow(parsed), 1)
  expect_equal(parsed$Cluster, "Cluster1")
  expect_equal(parsed$CellType, "T cell")
  expect_equal(parsed$Confidence, 0.95)
  expect_false(parsed$IsMixed)
})

test_that("marker input processing removes common low-value genes", {
  markers <- process_cell_data(list(Cluster1 = "CD3D, CD3E; MT-ATP6\nRPL13A CD8A"))$markers

  expect_equal(markers$Cluster1, c("CD3D", "CD3E", "CD8A"))
})

test_that("fallback ontology maps common aliases with provenance", {
  ontology <- create_fallback_ontology()
  mapping <- map_to_cell_ontology("NK cell", ontology)

  expect_equal(mapping$CL_ID, "CL:0000623")
  expect_true(mapping$MatchMethod %in% c("exact", "synonym"))
  expect_true(mapping$OntologyMatchScore > 0.9)
})

test_that("pancreas context disambiguates endocrine and ductal labels", {
  ontology <- create_fallback_ontology()

  expect_equal(
    map_to_cell_ontology("Alpha cell", ontology, tissue = "Pancreas")$CL_ID,
    "CL:0000171"
  )
  expect_equal(
    map_to_cell_ontology("Beta cell", ontology, tissue = "Pancreas")$CL_ID,
    "CL:0000169"
  )
  expect_equal(
    map_to_cell_ontology("Delta cell", ontology, tissue = "Pancreas")$CL_ID,
    "CL:0000173"
  )
  expect_equal(
    map_to_cell_ontology("Acinar cell", ontology, tissue = "Pancreas")$CL_ID,
    "CL:0002064"
  )
  expect_equal(
    map_to_cell_ontology("Ductal cell", ontology, tissue = "Pancreas")$CL_ID,
    "CL:0002079"
  )
})

test_that("brain context maps common neural and glial labels", {
  ontology <- create_fallback_ontology()

  expect_equal(
    map_to_cell_ontology("Neuron", ontology, tissue = "Brain")$CL_ID,
    "CL:2000029"
  )
  expect_equal(
    map_to_cell_ontology("Astrocytes", ontology, tissue = "Brain")$CL_ID,
    "CL:0000127"
  )
  expect_equal(
    map_to_cell_ontology("Microglia", ontology, tissue = "Brain")$CL_ID,
    "CL:0000129"
  )
  expect_equal(
    map_to_cell_ontology("OPC", ontology, tissue = "Brain")$CL_ID,
    "CL:0002453"
  )
  expect_equal(
    map_to_cell_ontology("Excitatory neuron", ontology, tissue = "Brain")$CL_ID,
    "CL:0000679"
  )
  expect_equal(
    map_to_cell_ontology("Inhibitory neuron", ontology, tissue = "Brain")$CL_ID,
    "CL:0000617"
  )
  expect_equal(
    map_to_cell_ontology("Pyramidal cell", ontology, tissue = "Brain")$CL_ID,
    "CL:0000598"
  )
  expect_equal(
    map_to_cell_ontology("Ependymal", ontology, tissue = "Brain")$CL_ID,
    "CL:0000065"
  )
})

test_that("lung context maps common epithelial, stromal, and immune labels", {
  ontology <- create_fallback_ontology()

  expect_equal(
    map_to_cell_ontology("AT1 cell", ontology, tissue = "Lung")$CL_ID,
    "CL:0002062"
  )
  expect_equal(
    map_to_cell_ontology("AT2 cell", ontology, tissue = "Lung")$CL_ID,
    "CL:0002063"
  )
  expect_equal(
    map_to_cell_ontology("Clara cell", ontology, tissue = "Lung")$CL_ID,
    "CL:0000158"
  )
  expect_equal(
    map_to_cell_ontology("Ciliated cell", ontology, tissue = "Lung")$CL_ID,
    "CL:1000271"
  )
  expect_equal(
    map_to_cell_ontology("Basal cell", ontology, tissue = "Lung")$CL_ID,
    "CL:0002633"
  )
  expect_equal(
    map_to_cell_ontology("Goblet cell", ontology, tissue = "Lung")$CL_ID,
    "CL:1000143"
  )
  expect_equal(
    map_to_cell_ontology("Endothelial cell", ontology, tissue = "Lung")$CL_ID,
    "CL:1001567"
  )
  expect_equal(
    map_to_cell_ontology("Fibroblast", ontology, tissue = "Lung")$CL_ID,
    "CL:0002553"
  )
  expect_equal(
    map_to_cell_ontology("Alveolar macrophage", ontology, tissue = "Lung")$CL_ID,
    "CL:0000583"
  )
  expect_equal(
    map_to_cell_ontology("RASC", ontology, tissue = "Lung")$CL_ID,
    "CL:4052031"
  )
  expect_equal(
    map_to_cell_ontology("Ionocyte", ontology, tissue = "Lung")$CL_ID,
    "CL:0017000"
  )
  expect_equal(
    map_to_cell_ontology("Smooth muscle cell", ontology, tissue = "Lung")$CL_ID,
    "CL:0019019"
  )
})

test_that("validation reports quality metrics", {
  validation <- validate_annotations(
    data.frame(
      Cluster = c("Cluster1", "Cluster2"),
      CellType = c("T cell", "Unknown"),
      Confidence = c(0.95, 0.4),
      CL_ID = c("CL:0000084", NA),
      IsMixed = c(FALSE, TRUE),
      stringsAsFactors = FALSE
    )
  )

  expect_named(validation, c(
    "valid", "issues", "warnings", "summary", "thresholds",
    "metadata", "timestamp", "quality_score"
  ))
  expect_equal(validation$summary$n_clusters, 2)
  expect_equal(validation$summary$ontology_coverage, 0.5)
})
