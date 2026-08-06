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

test_that("LLM JSON responses preserve ranked candidate labels", {
  response <- paste0(
    '{"annotations":[{"cluster":"Cluster1","cell_type":"T cell","confidence":0.9,',
    '"candidate_cell_types":["T cell","NK cell"],',
    '"is_mixed":false,"tissue_consistency":"expected","reasoning":"CD3 markers"}]}'
  )

  parsed <- parse_annotation_response(response)

  expect_equal(parsed$CandidateCellTypes, "T cell; NK cell")
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

test_that("evidence-adjusted confidence preserves strong marker-supported annotations", {
  annotations <- data.frame(
    Cluster = "Cluster1",
    CellType = "Beta cell",
    Confidence = 0.9,
    IsMixed = FALSE,
    TissueConsistency = "expected",
    CL_ID = "CL:0000169",
    OntologyLabel = "type B pancreatic cell",
    MatchMethod = "context_exact",
    OntologyMatchScore = 1,
    stringsAsFactors = FALSE
  )
  markers <- list(Cluster1 = c("INS", "IAPP", "MAFA", "PDX1"))

  calibrated <- calibrate_annotation_confidence(annotations, markers, tissue = "Pancreas")

  expect_equal(calibrated$LLMConfidence, 0.9)
  expect_equal(calibrated$EvidenceBestCellType, "beta cell")
  expect_false(calibrated$EvidenceConflict)
  expect_gte(calibrated$Confidence, 0.85)
  expect_equal(calibrated$ConfidenceMethod, "ontology_marker_calibrated")
})

test_that("evidence scoring flags marker conflicts for refinement", {
  annotations <- data.frame(
    Cluster = "Cluster1",
    CellType = "Macrophage",
    Confidence = 0.95,
    IsMixed = FALSE,
    TissueConsistency = "expected",
    CL_ID = "CL:0000235",
    OntologyLabel = "macrophage",
    MatchMethod = "exact",
    OntologyMatchScore = 1,
    Reasoning = "first-pass label",
    stringsAsFactors = FALSE
  )
  markers <- list(Cluster1 = c("INS", "IAPP", "MAFA", "PDX1"))

  scored <- score_annotation_evidence(annotations, markers, tissue = "Pancreas")
  prompt <- create_refinement_prompt(
    markers,
    cbind(annotations, scored[setdiff(names(scored), "Cluster")]),
    tissue = "Pancreas"
  )

  expect_equal(scored$EvidenceBestCellType, "beta cell")
  expect_true(scored$EvidenceConflict)
  expect_true(scored$RequiresRefinement)
  expect_match(prompt, "evidence-guided selective refinement", ignore.case = TRUE)
  expect_match(prompt, "INS", fixed = TRUE)
})

test_that("fixed-budget selector supports evidence, confidence, and full strategies", {
  annotations <- data.frame(
    Cluster = c("Cluster1", "Cluster2", "Cluster3"),
    CellType = c("Macrophage", "Beta cell", "T cell"),
    Confidence = c(0.95, 0.2, 0.7),
    IsMixed = FALSE,
    TissueConsistency = "expected",
    CL_ID = c("CL:0000235", "CL:0000169", "CL:0000084"),
    OntologyLabel = c("macrophage", "type B pancreatic cell", "T cell"),
    MatchMethod = c("exact", "context_exact", "exact"),
    OntologyMatchScore = 1,
    stringsAsFactors = FALSE
  )
  markers <- list(
    Cluster1 = c("INS", "IAPP", "MAFA", "PDX1"),
    Cluster2 = c("INS", "IAPP", "MAFA", "PDX1"),
    Cluster3 = c("CD3D", "CD3E", "TRAC")
  )

  scored <- calibrate_annotation_confidence(annotations, markers, tissue = "Pancreas")

  evidence_k <- select_refinement_candidates(scored, strategy = "evidence", budget = 1)
  confidence_k <- select_refinement_candidates(scored, strategy = "confidence", budget = 1)
  full <- select_refinement_candidates(scored, strategy = "full")

  expect_equal(evidence_k$Cluster, "Cluster1")
  expect_equal(evidence_k$SelectionStrategy, "evidence")
  expect_equal(confidence_k$Cluster, "Cluster2")
  expect_equal(full$Cluster, scored$Cluster)
  expect_equal(nrow(full), 3)
})

test_that("ontology-disabled evidence removes ontology-only conflict triggers", {
  annotations <- data.frame(
    Cluster = "Cluster1",
    CellType = "Beta cell",
    Confidence = 0.9,
    IsMixed = FALSE,
    TissueConsistency = "expected",
    CL_ID = "CL:low_quality_match",
    OntologyLabel = "poor ontology match",
    MatchMethod = "fuzzy",
    OntologyMatchScore = 0.1,
    stringsAsFactors = FALSE
  )
  markers <- list(Cluster1 = c("INS", "IAPP", "MAFA", "PDX1"))

  with_ontology <- score_annotation_evidence(
    annotations,
    markers,
    tissue = "Pancreas",
    use_ontology_evidence = TRUE
  )
  without_ontology <- score_annotation_evidence(
    annotations,
    markers,
    tissue = "Pancreas",
    use_ontology_evidence = FALSE
  )

  expect_true(with_ontology$EvidenceConflict)
  expect_true(with_ontology$RequiresRefinement)
  expect_equal(without_ontology$OntologyEvidenceScore, 0.5)
  expect_false(without_ontology$EvidenceConflict)
  expect_false(without_ontology$RequiresRefinement)
})

test_that("refinement provenance distinguishes reviewed rows from label changes", {
  first_pass <- data.frame(
    Cluster = c("Cluster1", "Cluster2"),
    CellType = c("Macrophage", "Beta cell"),
    Confidence = c(0.95, 0.88),
    stringsAsFactors = FALSE
  )
  refined <- data.frame(
    Cluster = c("Cluster1", "Cluster2"),
    CellType = c("Beta cell", "Beta cell"),
    Confidence = c(0.93, 0.90),
    stringsAsFactors = FALSE
  )

  summary <- .summarise_refinement_changes(first_pass, refined)
  final <- .add_refinement_provenance(
    annotations = refined,
    first_pass_annotations = first_pass,
    flagged_clusters = c("Cluster1", "Cluster2"),
    refined_clusters = c("Cluster1", "Cluster2"),
    label_changed_clusters = summary$label_changed_clusters,
    confidence_changed_clusters = summary$confidence_changed_clusters
  )

  expect_equal(summary$n_label_changed, 1)
  expect_equal(summary$n_confidence_changed, 2)
  expect_true(final$WasFlagged[1])
  expect_true(final$WasRefined[1])
  expect_true(final$RefinementChangedLabel[1])
  expect_false(final$RefinementChangedLabel[2])
  expect_equal(final$FirstPassCellType[1], "Macrophage")
})
