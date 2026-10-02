# Related work: single-cell foundation models, reverse translation and breast-cancer prognosis

*Written 2026-10-02 for the study in this repository: cell states found by scGPT, Geneformer and an
HVG-PCA baseline (highly variable genes, principal component analysis) in a breast-cancer single-cell
atlas, turned into marker-gene signatures, scored for overall survival (OS) in TCGA-BRCA (The Cancer
Genome Atlas breast cohort) and tested again in METABRIC (Molecular Taxonomy of Breast Cancer
International Consortium), with PAM50 (Prediction Analysis of Microarray, 50 genes) subtype as a
reference. Only
references whose metadata were checked against Crossref, arXiv or the publisher are listed; how each was
checked is at the end. Descriptions marked † go beyond the paper's abstract and rest on a full-text
reading made on 2026-09-30.*

## 1. Foundation models against simple baselines

**The two models tested here.** scGPT (Cui et al. 2024) is used as its whole-human checkpoint, which its
repository recommends for most applications and describes as pretrained on 33 million normal human
cells. Geneformer (Theodoris et al. 2023) is used as its V1 10M-parameter model, pretrained on
Genecorpus-30M, which excluded malignant cells and immortalised cell lines. Neither model saw a malignant
cell in pretraining, while the states that carry survival signal in this study are mostly malignant-cell
states.

**Benchmarks against simple baselines.** Zero-shot evaluation of the same two checkpoints found that simpler methods
can outperform them (Kedzierska et al. 2025); highly variable genes, scVI and Harmony beat both on
cell-type clustering†. L1-regularised logistic regression matched fine-tuned scBERT and scGPT on
cell-type annotation (Boiarsky et al. 2023; the peer-reviewed follow-up, Boiarsky et al. 2024, covers
scBERT only; Geneformer was not evaluated). For perturbation prediction, none of five foundation models
and two other deep models outperformed simple linear baselines (Ahlmann-Eltze et al. 2025); the
phrase "one PCA still rules them all" is the title of a perturbation benchmark (Bendidi et al. 2024); a
mean-of-training baseline outperformed scGPT and scFoundation (Csendes et al. 2025†); zero-shot
embeddings gave no consistent gain in PertEval-scFM (Wenteler et al. 2024†); and systematic variation
inflates the usual perturbation metrics (Viñas Torné et al. 2026). Among newer models, a parameter-free
linear pipeline matched TranscriptFormer on cell-type classification and exceeded it on
out-of-distribution tasks (Souza & Mehta 2026; it does not evaluate scGPT or Geneformer), and a
zero-shot benchmark of scGPT, SCimilarity, Universal Cell Embeddings (UCE) and TranscriptFormer found the baselines best for
annotation† (Gaballa et al. 2026). Broader benchmarks find that simpler models adapt more efficiently
to specific datasets (Wu J et al. 2025) and that foundation models do not consistently outperform
task-specific methods (Liu T et al. 2026). A pre-registered evaluation of Geneformer, scGPT and UCE with
capacity and shuffle controls is the closest methodological sibling of this study's control ladder
(Shoeibi & Yousefi 2026).

**Scaling studies.** UCE (Rosen et al. 2026), TranscriptFormer
(Pearce et al. 2026), scFoundation (Hao et al. 2024), CellFM, pretrained on 100 million human cells (Zeng
et al. 2025), and Cell2Sentence scaled to 27 billion parameters (Rizvi et al. 2025†) argue that larger
corpora and models learn transferable cell biology. Geneformer V2 was trained on about 104 million cells
(Chen H et al. 2026); its model card lists V2-316M as the current default and a cancer-tuned V2-104M
variant continually pretrained on about 14 million cancer transcriptomes. scGPT also releases a
pan-cancer checkpoint (5.7 million cells) and a continually pretrained checkpoint for zero-shot cell
embedding. A negative result here is a result about the light V1 Geneformer and whole-human scGPT
checkpoints; the cancer-adapted variants were not tested.

## 2. From single-cell states to bulk outcomes (reverse translation)

**Deconvolution and ecotypes.** CIBERSORT (Newman et al. 2015) and CIBERSORTx (Newman et al. 2019)
estimate cell-type abundance in bulk from single-cell references; BayesPrism (Chu et al. 2022) does so
with a Bayesian model that separates malignant-cell expression from the microenvironment. EcoTyper (Luca
et al. 2021) identified 69 cell states across 16 carcinoma types, most of them prognostic, and ten
multicellular communities. Dai et al. (2026) review the field; Li et al. (2026) benchmark deconvolution
methods on real bulk cohorts, including whether prognostic associations reproduce across cohorts. This
study scores signatures as the mean z-score of their genes, a simpler proxy that mixes state abundance
with tumour purity and proliferation.

**Phenotype-guided cell selection.** Scissor (Sun et al. 2022), scAB (Zhang et al. 2022), DEGAS (Johnson
et al. 2022), SCIPAC (Gan et al. 2024), scPAS (Xie et al. 2025) and scSurv (Mizukoshi et al. 2026) use the
bulk phenotype, survival included, to choose or weight single cells. Their associations are therefore
not independent tests of outcome. In this study state discovery never sees the outcome, which is what
makes an outcome-permutation null valid for the whole family of signatures.

**Breast cancer.** Wu et al. (2021) built a single-cell and spatial breast atlas, a single-cell intrinsic
subtype classifier (SCSubtype), and deconvolved large bulk cohorts into nine ecotypes with distinct
outcomes; 5,291 of this study's cells come from that dataset. Pang et al. (2024) found seven consensus
malignant-cell states and, by deconvolution, an immune-related state associated with better survival.
Chen A et al. (2026) integrated more than 600,000 cells from 138 patients and ran a survival analysis of
every subpopulation in TCGA, METABRIC and SCAN-B (Sweden Cancerome Analysis Network – Breast); the associations they report as holding across cohorts
are states that confer a survival advantage. Their CELLxGENE collection supplies four of this study's
eight datasets and 31,552 of its 50,002 cells (63%). Malignant cells cluster by patient in single-cell
data (Tirosh et al. 2016), and recurrent malignant programmes are shared across tumours (Gavish et al.
2023; Barkley et al. 2022); several states in this study are single-donor malignant programmes.

## 3. Random signatures, subtype and proliferation

More than 90% of random signatures with over 100 genes predict breast-cancer outcome, 60% of 47 published
signatures were no better than random signatures of the same size, and adjusting for a proliferation
metagene removed almost all of these associations (Venet et al. 2011). Random gene sets predict survival
in most TCGA cancer types, and there a proliferation signature does not remove the effect while
data-driven sub-classification reduces it (Shimoni 2018). Significance Analysis of Prognostic Signatures
(SAPS) requires a gene set to beat random gene sets and to be enriched for individually prognostic genes
(Beck et al. 2013). Published breast signatures collapse onto proliferation, oestrogen-receptor and HER2
modules, with proliferation carrying most of the prognostic information (Wirapati et al. 2008; Sotiriou
et al. 2006); a single proliferation gene performs about as well as complex models (Haibe-Kains et al.
2008); signatures with little gene overlap give concordant predictions (Fan et al. 2006); and the ease of
finding prognostic signatures follows from confounding with subtype and clinical variables (Tofigh et al.
2014†). Which processes are prognostic depends on subtype: proliferation in oestrogen-receptor-positive, HER2-negative
(ER+/HER2−) disease, immune response in ER−/HER2− disease (Desmedt et al. 2008). Immune and differentiated states are often
protective (Gentles et al. 2015; Ali et al. 2016; Luca et al. 2021), and tumour purity confounds bulk
analyses (Aran et al. 2015).

How the study uses this: the floor is Venet's random-signature control, matched on size and on mean
expression bin as in Tirosh et al. (2016); the permutation null and the floor play the roles of SAPS's
tests against outcome permutation and against random gene sets, with a family-wise correction for the
best of several hundred signatures. The 11-gene proliferation score of PAM50 (Nielsen et al. 2010;
Parker et al. 2009) is the published reference signature. P3, prognosis within PAM50 subtype, is the
control for rediscovered subtype. The design read the C-index in the risk direction only, although the
literature above expects protective states; the post-hoc sensitivity analysis (DESIGN §11) reads each
signature in its own direction.

## 4. Foundation models and patient outcomes

The closest benchmark evaluated twelve single-cell foundation models and three baselines on seven cancer
tasks, including subtype classification and treatment response, and found limited advantages over the
simpler baselines for clinical and biological outcomes (Roman et al. 2025); its patients are classified
within single-cell cohorts, and it has no survival task†. Patient-level representations learned from
single-cell cohorts (Liu T et al. 2026, PaSCient) and drug-response benchmarks (Wang Q et al. 2025) take
the same within-cohort route.

A second route embeds the bulk tumour profile itself. scFoundation embeddings of TCGA bulk samples across
25 cancer types, combined with gene expression and clinical variables, reached a mean C-index of 0.724
and outperformed single-modality models (Liu W et al. 2026). CancerFoundation, trained on malignant cells,
proposes bulk survival prediction as a downstream task (Theus et al. 2024); Geneformer has been
fine-tuned on TCGA bulk profiles for overall survival (Mellors et al. 2024†); GeneBag fine-tunes a cell
model on bulk data (Liang et al. 2024); BulkRNABert (Gélard et al. 2025) and BulkFormer (Kang et al. 2026)
are pretrained on bulk transcriptomes. On TCGA survival, baseline methods matched or exceeded deep
representation-learning methods (Gross et al. 2024). A foundation-model-assisted breast signature has been
built and evaluated without a random-gene control (Liu W, Wu W et al. 2026†).

## 5. What is new here, and what is not

Not new: carrying single-cell states into bulk cohorts to test prognosis (Wu et al. 2021; Luca et al.
2021; Chen A et al. 2026); random-signature and proliferation controls (Venet et al. 2011; Beck et al.
2013); testing foundation models against simple baselines on clinical tasks (Roman et al. 2025; Gross et
al. 2024; Liu W et al. 2026).

We found no prior study that combines the following: a foundation model used only to define cell states
in a tumour atlas; those states carried as interpretable marker signatures into a disjoint bulk cohort
with survival follow-up; a linear HVG-PCA arm run through the identical clustering, signature and
scoring pipeline; every C-index read against a family-wise permutation null and a random-gene floor
matched on size and expression; predictions committed before the data; and a replication in an
independent cohort with signatures frozen before its outcomes were read. The contribution is
methodological and the answer agrees with the benchmarks above: in TCGA no representation's best
signature clears its family-wise null, and in METABRIC every foundation-model-minus-baseline margin is
negative. The protective luminal states that HVG-PCA and scGPT find in the post-hoc reading are in line
with Chen A et al. (2026), whose cross-cohort associations are protective.

## 6. Data sources and software

TCGA-BRCA (Cancer Genome Atlas Network 2012) through UCSC Xena (Goldman et al. 2020), with survival
endpoints from the TCGA Pan-Cancer Clinical Data Resource (Liu J et al. 2018), which flags breast-cancer
OS as needing longer follow-up†. METABRIC (Curtis et al. 2012; Pereira et al. 2016) from the cBioPortal
datahub. The atlas comes from CZ CELLxGENE Discover (CZI Cell Science Program et al. 2025); its eight
datasets were published by Chen A et al. (2026; four datasets), Wu et al. (2021), Gondal et al. (2025;
integrated immune-checkpoint-blockade data whose breast donors come from Bassez et al. 2021), Guimarães et
al. (2024; a multi-tissue atlas whose breast donors come from Qian et al. 2020) and Klughammer et al.
(2024; metastatic biopsies). Concordance index: Harrell et al. (1982). Clustering: Leiden (Traag et al.
2019) in Scanpy (Wolf et al. 2018). Survival models: lifelines (Davidson-Pilon 2019).

## 7. Corrections to the earlier version of this note and of the project page

- Kedzierska et al. 2025 was cited with its preprint title; the Genome Biology title is used now.
- Boiarsky et al. 2023: the last author is D. Sontag; the paper compares logistic regression with
  fine-tuned scBERT and scGPT, not Geneformer, on four datasets.
- "One PCA still rules them all" was attributed to Kedzierska et al.; it is the title of Bendidi et al.
  2024.
- Souza & Mehta 2026 compare a linear pipeline with TranscriptFormer; they do not evaluate scGPT,
  Geneformer or UCE, or perturbation tasks.
- Rosen et al. 2026: the third author is A. Agrawal.
- CZ CELLxGENE Discover: the corporate author is the CZI Cell Science Program.
- The earlier note described a different design (layer-resolved probing with UCE and scVI arms); it has
  been replaced.

## References

1. Ahlmann-Eltze C, Huber W, Anders S. Deep-learning-based gene perturbation effect prediction does not yet outperform simple linear baselines. *Nature Methods* 22:1657–1661 (2025). doi:10.1038/s41592-025-02772-6
2. Ali HR, Chlon L, Pharoah PDP, Markowetz F, Caldas C. Patterns of immune infiltration in breast cancer and their clinical implications: a gene-expression-based retrospective study. *PLoS Medicine* 13:e1002194 (2016). doi:10.1371/journal.pmed.1002194
3. Aran D, Sirota M, Butte AJ. Systematic pan-cancer analysis of tumour purity. *Nature Communications* 6:8971 (2015). doi:10.1038/ncomms9971
4. Barkley D, Moncada R, Pour M, et al. Cancer cell states recur across tumor types and form specific interactions with the tumor microenvironment. *Nature Genetics* 54:1192–1201 (2022). doi:10.1038/s41588-022-01141-9
5. Bassez A, Vos H, Van Dyck L, et al. A single-cell map of intratumoral changes during anti-PD1 treatment of patients with breast cancer. *Nature Medicine* 27:820–832 (2021). doi:10.1038/s41591-021-01323-8
6. Beck AH, Knoblauch NW, Hefti MM, et al. Significance Analysis of Prognostic Signatures. *PLoS Computational Biology* 9:e1002875 (2013). doi:10.1371/journal.pcbi.1002875
7. Bendidi I, Whitfield S, Kenyon-Dean K, Ben Yedder H, El Mesbahi Y, Noutahi E, Denton AK. Benchmarking transcriptomics foundation models for perturbation analysis: one PCA still rules them all. arXiv:2410.13956 (2024).
8. Boiarsky R, Singh N, Buendia A, Getz G, Sontag D. A deep dive into single-cell RNA sequencing foundation models. bioRxiv (2023). doi:10.1101/2023.10.19.563100
9. Boiarsky R, Singh NM, Buendia A, Amini AP, Getz G, Sontag D. Deeper evaluation of a single-cell foundation model. *Nature Machine Intelligence* 6:1443–1446 (2024). doi:10.1038/s42256-024-00949-w
10. Cancer Genome Atlas Network. Comprehensive molecular portraits of human breast tumours. *Nature* 490:61–70 (2012). doi:10.1038/nature11412
11. Chen A, Kroehling L, Ennis CS, Denis GV, Monti S. A highly resolved integrated single-cell atlas of human breast cancers. *NAR Genomics and Bioinformatics* 8:lqaf217 (2026). doi:10.1093/nargab/lqaf217
12. Chen H, Venkatesh MS, Gómez Ortega J, Mahesh SV, Nandi TN, Madduri RK, Pelka K, Theodoris CV. Scaling and quantization of large-scale foundation model enables resource-efficient predictions in network biology. *Nature Computational Science* 6:450–463 (2026). doi:10.1038/s43588-026-00972-4
13. Chu T, Wang Z, Pe'er D, Danko CG. Cell type and gene expression deconvolution with BayesPrism enables Bayesian integrative analysis across bulk and single-cell RNA sequencing in oncology. *Nature Cancer* 3:505–517 (2022). doi:10.1038/s43018-022-00356-3
14. Csendes G, Sanz G, Szalay KZ, Szalai B. Benchmarking foundation cell models for post-perturbation RNA-seq prediction. *BMC Genomics* 26:393 (2025). doi:10.1186/s12864-025-11600-2
15. Cui H, Wang C, Maan H, Pang K, Luo F, Duan N, Wang B. scGPT: toward building a foundation model for single-cell multi-omics using generative AI. *Nature Methods* 21:1470–1480 (2024). doi:10.1038/s41592-024-02201-0
16. Curtis C, Shah SP, Chin SF, et al. The genomic and transcriptomic architecture of 2,000 breast tumours reveals novel subgroups. *Nature* 486:346–352 (2012). doi:10.1038/nature10983
17. CZI Cell Science Program, Abdulla S, Aevermann B, et al. CZ CELLxGENE Discover: a single-cell data platform for scalable exploration, analysis and modeling of aggregated data. *Nucleic Acids Research* 53:D886–D900 (2025). doi:10.1093/nar/gkae1142
18. Dai Y, Guo S, Pan Y, Castignani C, Montierth MD, Van Loo P, Wang W. A guide to transcriptomic deconvolution in cancer. *Nature Reviews Cancer* 26:84–103 (2026). doi:10.1038/s41568-025-00886-9
19. Davidson-Pilon C. lifelines: survival analysis in Python. *Journal of Open Source Software* 4:1317 (2019). doi:10.21105/joss.01317
20. Desmedt C, Haibe-Kains B, Wirapati P, et al. Biological processes associated with breast cancer clinical outcome depend on the molecular subtypes. *Clinical Cancer Research* 14:5158–5165 (2008). doi:10.1158/1078-0432.CCR-07-4756
21. Fan C, Oh DS, Wessels L, Weigelt B, Nuyten DSA, Nobel AB, van't Veer LJ, Perou CM. Concordance among gene-expression-based predictors for breast cancer. *New England Journal of Medicine* 355:560–569 (2006). doi:10.1056/NEJMoa052933
22. Gaballa Y, Ahmed S, Abdelaal T. Benchmarking single-cell foundation models in a zero-shot setting. bioRxiv (2026). doi:10.64898/2026.08.03.739553
23. Gan D, Zhu Y, Lu X, Li J. SCIPAC: quantitative estimation of cell-phenotype associations. *Genome Biology* 25:119 (2024). doi:10.1186/s13059-024-03263-1
24. Gavish A, Tyler M, Greenwald AC, et al. Hallmarks of transcriptional intratumour heterogeneity across a thousand tumours. *Nature* 618:598–606 (2023). doi:10.1038/s41586-023-06130-4
25. Gélard M, Richard G, Pierrot T, Cournède P-H. BulkRNABert: cancer prognosis from bulk RNA-seq based language models. *Proceedings of the 4th Machine Learning for Health Symposium*, PMLR 259:384–400 (2025). bioRxiv doi:10.1101/2024.06.18.599483
26. Gentles AJ, Newman AM, Liu CL, et al. The prognostic landscape of genes and infiltrating immune cells across human cancers. *Nature Medicine* 21:938–945 (2015). doi:10.1038/nm.3909
27. Goldman MJ, Craft B, Hastie M, et al. Visualizing and interpreting cancer genomics data via the Xena platform. *Nature Biotechnology* 38:675–678 (2020). doi:10.1038/s41587-020-0546-8
28. Gondal M, Cieslik M, Chinnaiyan A. Integrated cancer cell-specific single-cell RNA-seq datasets of immune checkpoint blockade-treated patients. *Scientific Data* 12:139 (2025). doi:10.1038/s41597-025-04381-6
29. Gross B, Dauvin A, Cabeli V, Kmetzsch V, El Khoury J, Dissez G, et al. Robust evaluation of deep learning-based representation methods for survival and gene essentiality prediction on bulk RNA-seq data. *Scientific Reports* 14:17064 (2024). doi:10.1038/s41598-024-67023-8
30. Guimarães G, Maklouf G, Teixeira C, et al. Single-cell resolution characterization of myeloid-derived cell states with implication in cancer outcome. *Nature Communications* 15:5694 (2024). doi:10.1038/s41467-024-49916-4
31. Haibe-Kains B, Desmedt C, Sotiriou C, Bontempi G. A comparative study of survival models for breast cancer prognostication based on microarray data: does a single gene beat them all? *Bioinformatics* 24:2200–2208 (2008). doi:10.1093/bioinformatics/btn374
32. Hao M, Gong J, Zeng X, et al. Large-scale foundation model on single-cell transcriptomics. *Nature Methods* 21:1481–1491 (2024). doi:10.1038/s41592-024-02305-7
33. Harrell FE, Califf RM, Pryor DB, Lee KL, Rosati RA. Evaluating the yield of medical tests. *JAMA* 247:2543–2546 (1982). doi:10.1001/jama.1982.03320430047030
34. Johnson TS, Yu CY, Huang Z, et al. Diagnostic Evidence GAuge of Single cells (DEGAS): a flexible deep transfer learning framework for prioritizing cells in relation to disease. *Genome Medicine* 14:11 (2022). doi:10.1186/s13073-022-01012-2
35. Kang B, Fan R, Yi M, Cui C, Cui Q. BulkFormer: a large-scale foundation model for bulk transcriptomes. *Cell Systems* 17:101657 (2026). doi:10.1016/j.cels.2026.101657
36. Kedzierska KZ, Crawford L, Amini AP, Lu AX. Zero-shot evaluation reveals limitations of single-cell foundation models. *Genome Biology* 26:101 (2025). doi:10.1186/s13059-025-03574-x
37. Klughammer J, Abravanel DL, Segerstolpe Å, et al. A multi-modal single-cell and spatial expression map of metastatic breast cancer biopsies across clinicopathological features. *Nature Medicine* 30:3236–3249 (2024). doi:10.1038/s41591-024-03215-z
38. Li M, Su Y, Tang Y, Lee Y, Tian W. Evaluating deconvolution methods using real bulk RNA-expression data for robust prognostic insights across cancer types. *Genome Biology* 27:38 (2026). doi:10.1186/s13059-026-03942-1
39. Liang Y, Li D, Xu AG, Shao Y, Tang K. GeneBag: training a cell foundation model for broad-spectrum cancer diagnosis and prognosis with bulk RNA-seq data. bioRxiv (2024). doi:10.1101/2024.06.27.601098
40. Liu J, Lichtenberg T, Hoadley KA, et al. An integrated TCGA pan-cancer clinical data resource to drive high-quality survival outcome analytics. *Cell* 173:400–416.e11 (2018). doi:10.1016/j.cell.2018.02.052
41. Liu T, De Brouwer E, Verma A, et al. Learning multi-cellular representations of single-cell transcriptomics data enables characterization of patient-level disease states. *Cell Systems* 17:101570 (2026). doi:10.1016/j.cels.2026.101570
42. Liu T, Li K, Wang Y, Li H, Zhao H. Evaluating the utilities of foundation models in single-cell data analysis. *Advanced Science* 13:e14490 (2026). doi:10.1002/advs.202514490
43. Liu W, Wang Q, Long L, Wang W. Leveraging single-cell foundation models for accurate survival outcome prediction. *Bioinformatics Advances* 6:vbag076 (2026). doi:10.1093/bioadv/vbag076
44. Liu W, Wu W, Chen S, Chen K, Li X. AI-assisted multimodal transcriptomic analysis identifies a senescence-related prognostic signature and characterizes ADGRF5-associated malignant phenotypes in breast cancer. *Cells* 15:1589 (2026). doi:10.3390/cells15171589
45. Luca BA, Steen CB, Matusiak M, et al. Atlas of clinically distinct cell states and ecosystems across human solid tumors. *Cell* 184:5482–5496.e28 (2021). doi:10.1016/j.cell.2021.09.014
46. Mellors T, Schneider M, Spitmann M. A transformer-based approach to survival outcome prediction. bioRxiv (2024). doi:10.1101/2024.11.03.621674
47. Mizukoshi C, Kojima Y, Hayashi S, Abe K, Kasugai D, Shimamura T. scSurv: a deep generative model for single-cell survival analysis. *Bioinformatics* 42:btaf671 (2026). doi:10.1093/bioinformatics/btaf671
48. Newman AM, Liu CL, Green MR, et al. Robust enumeration of cell subsets from tissue expression profiles. *Nature Methods* 12:453–457 (2015). doi:10.1038/nmeth.3337
49. Newman AM, Steen CB, Liu CL, et al. Determining cell type abundance and expression from bulk tissues with digital cytometry. *Nature Biotechnology* 37:773–782 (2019). doi:10.1038/s41587-019-0114-2
50. Nielsen TO, Parker JS, Leung S, et al. A comparison of PAM50 intrinsic subtyping with immunohistochemistry and clinical prognostic factors in tamoxifen-treated estrogen receptor-positive breast cancer. *Clinical Cancer Research* 16:5222–5232 (2010). doi:10.1158/1078-0432.CCR-10-1282
51. Pang L, Xiang F, Yang H, et al. Single-cell integrative analysis reveals consensus cancer cell states and clinical relevance in breast cancer. *Scientific Data* 11:289 (2024). doi:10.1038/s41597-024-03127-0
52. Parker JS, Mullins M, Cheang MCU, et al. Supervised risk predictor of breast cancer based on intrinsic subtypes. *Journal of Clinical Oncology* 27:1160–1167 (2009). doi:10.1200/JCO.2008.18.1370
53. Pearce JD, Simmonds SE, Mahmoudabadi G, et al. TranscriptFormer: a generative cell atlas across 1.5 billion years of evolution. *Science* 393:aec8514 (2026). doi:10.1126/science.aec8514
54. Pereira B, Chin SF, Rueda OM, et al. The somatic mutation profiles of 2,433 breast cancers refine their genomic and transcriptomic landscapes. *Nature Communications* 7:11479 (2016). doi:10.1038/ncomms11479
55. Qian J, Olbrecht S, Boeckx B, et al. A pan-cancer blueprint of the heterogeneous tumor microenvironment revealed by single-cell profiling. *Cell Research* 30:745–762 (2020). doi:10.1038/s41422-020-0355-0
56. Rizvi SA, Levine D, Patel A, et al. Scaling large language models for next-generation single-cell analysis. bioRxiv (2025). doi:10.1101/2025.04.14.648850
57. Roman A, Johri S, Conci R, Van Allen EM, Elmarakeby H. Empirical evaluation of single-cell foundation models for predicting cancer outcomes. bioRxiv (2025; v2 2026). doi:10.1101/2025.10.31.685892
58. Rosen Y, Roohani Y, Agrawal A, Samotorčan L, Tabula Sapiens Consortium, Quake SR, Leskovec J. Universal cell embedding provides a foundation model for cell biology. *Nature* 656:183–191 (2026). doi:10.1038/s41586-026-10689-z
59. Shimoni Y. Association between expression of random gene sets and survival is evident in multiple cancer types and may be explained by sub-classification. *PLoS Computational Biology* 14:e1006026 (2018). doi:10.1371/journal.pcbi.1006026
60. Shoeibi M, Yousefi N. Pre-registered external evaluation yields a consistent partial-replication category across three transcriptomic foundation models. arXiv:2608.26170 (2026).
61. Sotiriou C, Wirapati P, Loi S, et al. Gene expression profiling in breast cancer: understanding the molecular basis of histologic grade to improve prognosis. *Journal of the National Cancer Institute* 98:262–272 (2006). doi:10.1093/jnci/djj052
62. Souza H, Mehta P. Parameter-free representations outperform single-cell foundation models on downstream benchmarks. arXiv:2602.16696 (2026).
63. Sun D, Guan X, Moran AE, et al. Identifying phenotype-associated subpopulations by integrating bulk and single-cell sequencing data. *Nature Biotechnology* 40:527–538 (2022). doi:10.1038/s41587-021-01091-3
64. Theodoris CV, Xiao L, Chopra A, et al. Transfer learning enables predictions in network biology. *Nature* 618:616–624 (2023). doi:10.1038/s41586-023-06139-9
65. Theus A, Barkmann F, Wissel D, Boeva V. CancerFoundation: a single-cell RNA sequencing foundation model to decipher drug resistance in cancer. bioRxiv (2024). doi:10.1101/2024.11.01.621087
66. Tirosh I, Izar B, Prakadan SM, et al. Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq. *Science* 352:189–196 (2016). doi:10.1126/science.aad0501
67. Tofigh A, Suderman M, Paquet ER, et al. The prognostic ease and difficulty of invasive breast carcinoma. *Cell Reports* 9:129–142 (2014). doi:10.1016/j.celrep.2014.08.073
68. Traag VA, Waltman L, van Eck NJ. From Louvain to Leiden: guaranteeing well-connected communities. *Scientific Reports* 9:5233 (2019). doi:10.1038/s41598-019-41695-z
69. Venet D, Dumont JE, Detours V. Most random gene expression signatures are significantly associated with breast cancer outcome. *PLoS Computational Biology* 7:e1002240 (2011). doi:10.1371/journal.pcbi.1002240
70. Viñas Torné R, Wiatrak M, Piran Z, Fan S, Jiang L, Teichmann SA, Nitzan M, Brbić M. Systema: a framework for evaluating genetic perturbation response prediction beyond systematic variation. *Nature Biotechnology* 44:1050–1059 (2026). doi:10.1038/s41587-025-02777-8
71. Wang Q, Pan Y, Zhou M, Tang Z, Wang Y, Wang G, Song Q. scDrugMap: benchmarking large foundation models for drug response prediction. *Nature Communications* 17:730 (2026). doi:10.1038/s41467-025-67481-2
72. Wenteler A, Occhetta M, Branson N, et al. PertEval-scFM: benchmarking single-cell foundation models for perturbation effect prediction. bioRxiv (2024). doi:10.1101/2024.10.02.616248
73. Wirapati P, Sotiriou C, Kunkel S, et al. Meta-analysis of gene expression profiles in breast cancer: toward a unified understanding of breast cancer subtyping and prognosis signatures. *Breast Cancer Research* 10:R65 (2008). doi:10.1186/bcr2124
74. Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology* 19:15 (2018). doi:10.1186/s13059-017-1382-0
75. Wu J, Ye Q, Wang Y, et al. Biology-driven insights into the power of single-cell foundation models. *Genome Biology* 26:334 (2025). doi:10.1186/s13059-025-03781-6
76. Wu SZ, Al-Eryani G, Roden DL, et al. A single-cell and spatially resolved atlas of human breast cancers. *Nature Genetics* 53:1334–1347 (2021). doi:10.1038/s41588-021-00911-1
77. Xie A, Wang H, Zhao J, Wang Z, Xu J, Xu Y. scPAS: single-cell phenotype-associated subpopulation identifier. *Briefings in Bioinformatics* 26:bbae655 (2025). doi:10.1093/bib/bbae655
78. Zeng Y, Xie J, Shangguan N, et al. CellFM: a large-scale foundation model pre-trained on transcriptomics of 100 million human cells. *Nature Communications* 16:4679 (2025). doi:10.1038/s41467-025-59926-5
79. Zhang Q, Jin S, Zou X. scAB detects multiresolution cell states with clinical significance by integrating single-cell genomics and bulk sequencing data. *Nucleic Acids Research* 50:12112–12130 (2022). doi:10.1093/nar/gkac1109

## How the references were verified

- **Metadata.** Every DOI above was resolved through the Crossref API on 2026-10-02 (76 DOIs, all
  found); title, journal, volume, pages or article number, year, first and last author were compared with
  the entry. The three arXiv preprints (Bendidi et al. 2024, Souza & Mehta 2026, Shoeibi & Yousefi 2026)
  were checked against their arXiv abstract pages. Where a paper appeared online one year and in an issue
  the next, the issue year is given.
- **Content.** Statements about what a paper found were checked against its abstract (Europe PMC) for
  Venet et al. 2011, Shimoni 2018, Beck et al. 2013, Wirapati et al. 2008, Desmedt et al. 2008, Kedzierska
  et al. 2025, Ahlmann-Eltze et al. 2025, Luca et al. 2021, Wu et al. 2021, Pang et al. 2024, Sun et al.
  2022, Chen A et al. 2026, Roman et al. 2025, Liu W et al. 2026 and Gross et al. 2024. Statements marked
  † rest on a full-text reading on 2026-09-30.
- **Models and data.** The scGPT checkpoint descriptions are from the scGPT repository README; the
  Geneformer corpus and checkpoint descriptions from the Geneformer model card and the Genecorpus-30M
  dataset card; the mapping of atlas datasets to publications from the CELLxGENE curation API
  (collection and publication DOI per dataset), checked on 2026-10-02.
- **Additions beyond the earlier reference audit.** Gondal et al. 2025 and Guimarães et al. 2024 were not in
  the audit's lists, which credit the same atlas datasets to their source studies (Bassez et al. 2021,
  Qian et al. 2020). They are cited because the CELLxGENE curation API names them as the publications of
  datasets `7b20c613` ("Integrated cancer cell-specific single-cell RNA-seq datasets of immune checkpoint
  blockade-treated patients") and `b617ee1b` ("A multi-tissue single-cell tumor microenvironment atlas");
  their metadata were checked against Crossref on 2026-10-02 like every other entry.
- Software repositories without a paper and references that could not be checked are not listed.
