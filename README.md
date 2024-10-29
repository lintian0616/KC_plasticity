Alternations in inflammatory macrophage niche drive phenotypic and functional plasticity of Kupffer cells
=========

This repository contains the source data and source codes for bioinformatic analysis. The paper was published in [Nature Communications](https://www.nature.com/articles/s41467-024-53659-7).

This study explores how liver-resident macrophages, known as Kupffer cells, change in the presence of liver metastatic nodules. The research demonstrates that Kupffer cells can be reprogrammed into an inflammatory state due to an inflammatory environment, which displaces them with monocyte-derived macrophages through a process called epigenetic reprogramming. This transformation is thought to contribute to the immune-suppressive environment that hinders liver metastasis treatment.

| Directory  | File  | Description |
|:------------- |:---------------|:-------------|
| CITEseq      | CITEseq_FeaturePlot.Rmd | related to **Fig. 1** |
| WPRE_scRNAseq| WPRE_DataVisualization.Rmd | related to **Fig. 4** |
| SLAM_ITseq | run\_full\_pipeline\_v2.sh & SLAMseq_Manuscript.Rmd |            related to **Supplementary Fig. 6** |
| EpigeneticReprogramming | CUT\_Tag\_Pipeline.md | related to **Fig. 5** |
| scTCR-seq | three R script files | related to **Supplementary Fig. 10** |

## Source Code/Data Download

You can use `git` command tool to clone the whole repository.

```
git clone https://github.com/lintian0616/KC_plasticity.git
```

Alternatively, you can click **Download ZIP** and later unzip the file in your local computer.

## Processed Data Download

The processed RData objects have been uploaded to [Zendo #10809097](https://zenodo.org/records/10809097) (macrophage scRNA-seq datasets) and [Zendo #12794202](https://zenodo.org/records/12794202) (CD8+ T cell scTCR-seq dataset).

Should you have any question, please do not hesitate to contact [Lin Tian](https://www.tianlab.info/)  (tianlin@sysucc.org.cn).