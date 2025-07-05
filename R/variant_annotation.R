#' Extract nucleotide change from variant string
#'
#' @description
#' Given a variant in the format "<pos>_<ref>_<alt>", returns the substring representing
#' the reference-to-alternate change (e.g. "A_T").
#'
#' @param variant Character vector of variant strings ("<pos>_<ref>_<alt>").
#' @return Character vector of nucleotide changes ("<ref>_<alt>").
#' @examples
#' add_changes("150_A_G")  # returns "A_G"
#' @export
add_changes <-function(variant){
    changes<-sub("^\\d+_", "",variant)
    return(changes)
}

#' Classify nucleotide change as transition or transversion
#'
#' @description
#' Determines whether each nucleotide change is a transition (purine-to-purine
#' or pyrimidine-to-pyrimidine) or a transversion (purine-to-pyrimidine or vice versa).
#'
#' @param changes Character vector of nucleotide changes ("<ref>_<alt>").
#' @return Character vector with values "transition" or "transversion".
#' @examples
#' add_types(c("C_T", "A_C"))  # returns c("transition", "transversion")
#' @export
add_types <- function(changes){
    types<- ifelse(changes %in% c("C_T","G_A","T_C","A_G"),"transition","transversion")
    return(types)
}

#' Annotate variant annotation with base-change and type
#'
#' @description
#' Adds two new columns to a variant annotation table:
#' 1. \code{changes}: the nucleotide substitution (e.g. "C_T").
#' 2. \code{types}: whether the substitution is a transition or transversion.
#'
#' @param variant_annotation A data.frame or tibble containing a column
#'   named \code{Variants} ("<pos>_<ref>_<alt>").
#' @return The input \code{variant_annotation} with added \code{changes} and \code{types} columns.
#' @examples
#' df <- data.frame(Variants = c("100_C_T", "200_A_C"))
#' Annotate_base_change(df)
#' @importFrom dplyr mutate
#' @export
Annotate_base_change <- function(variant_annotation){
    variant_annotation<- variant_annotation %>% mutate(changes=add_changes(Variants)) %>% mutate(types=add_types(changes))
    return(variant_annotation)
}


#' Annotate variant table with population frequency and haplogroup statistics
#'
#' @description
#' Adds mitochondrial variant information from \code{mitomap_freq} and 
#' haplogroup marker counts from \code{haplogroup_markers} to a variant annotation
#' data frame.  data. It drops any existing columns that would collide, then merges in:
#' \itemize{
#'   \item{\code{Locus}}{Gene name, tRNA, or other feature from \code{mitomap_freq}.}
#'   \item{\code{RSRS50}}{Logical flag indicating whether the variant is an ancestral (RSRS50) variant.}
#'   \item{\code{freq_all, freq_african, freq_asian, freq_eurasian}}{Mitochondrial variant frequencies in the overall, African, Asian, and Eurasian populations, respectively.}
#'   \item{\code{n_haplos}}{Number of population‐level ancestry branches (haplogroups) that consider this variant a marker (posterior >80\%).}
#' }
#'
#' @param variant_annotation A data.frame or tibble containing at least columns \code{Variants} and \code{CellN}.
#'
#' @return The input \code{variant_annotation}, augmented with:
#' \describe{
#'   \item{\code{Locus}}{Feature name (gene/tRNA/etc).}
#'   \item{\code{RSRS50}}{Ancestral variant flag (logical).}
#'   \item{\code{freq_all}}{Overall variant frequency.}
#'   \item{\code{freq_african}}{Variant frequency in African population.}
#'   \item{\code{freq_asian}}{Variant frequency in Asian population.}
#'   \item{\code{freq_eurasian}}{Variant frequency in Eurasian population.}
#'   \item{\code{n_haplos}}{Count of haplogroups (ancestry branches) marking this variant (>80\% posterior).}
#' }
#'
#' @importFrom dplyr select one_of left_join arrange desc group_by summarize
#' @export
annotate_variants_population_stats <- function(variant_annotation){ 
    data(mitomap_freq)
    data(haplogroup_markers)
    mitomap_freq <-mitomap_freq[,c("ID","Locus", "RSRS50","Overall.Variant.Frequency.in.Sequence.Set", "Frequency.in.Lineage.L", "Frequency.in.Lineage.M", "Frequency.in.Lineage.N")]
    names(mitomap_freq)[4:7] <-c("freq_all","freq_african","freq_asian","freq_eurasian")
    haplogroup_markers_summary <- haplogroup_markers %>%
      group_by(ID) %>%
      summarize(
        haplo_info = paste0(
          paste(haplo, collapse = "_"),  # collapse haplo by “_”
          ";",
          paste(freq, collapse = "_")    # collapse freq by “_”
        ),
        .groups = "drop",
        n_haplos = n()
      ) 
    variant_annotation <- variant_annotation %>% 
                        select(-one_of(unique(c(names(mitomap_freq),names(haplogroup_markers_summary))))) %>%
                        left_join(mitomap_freq, by = c("Variants" = "ID")) %>% 
                        left_join(haplogroup_markers_summary[, c("ID","n_haplos")], by = c("Variants" = "ID"))

    return(variant_annotation)
}


#' Annotate variants with mitochondrial blacklist region flag
#'
#' @description
#' Flags variants whose genomic position falls within known problematic
#' regions of the mitochondrial genome (blacklist regions).
#'
#' @details
#' The blacklist intervals correspond to regions in the revised Cambridge Reference
#' Sequence (rCRS) prone to misalignment:
#' \enumerate{
#'   \item Misalignment due to \code{ACCCCCC T C C C C C C} (rCRS 302–315)
#'   \item Misalignment due to \code{G C A C A C A C A C A C C} (rCRS 513–525)
#'   \item Misalignment due to \code{3107N} ambiguity in rCRS (rCRS 3105–3109)
#'   \item Misalignment due to \code{ACCCCC} homopolymer (rCRS 16182–16187)
#' }
#'
#' @param variant_annotation A data.frame or tibble containing at least a column
#'   named \code{Variants}, where each entry is of the form \code{"<pos>_<ref>_<alt>"}.
#'
#' @return The input \code{variant_annotation} with an added column:
#' \describe{
#'   \item{\code{blacklist_region}}{Character; \code{"blacklist_region"} if the variant
#'     position is in one of the predefined blacklist intervals, otherwise \code{""}.}
#' }
#'
#' @export
annotate_variants_blacklist <- function(variant_annotation){
    black_list_regions<-c(302:315, 513:525, 3105:3109, 16182:16187)
    pos<- as.numeric(strsplit(variant_annotation$Variants, "_") %>% sapply(function(x){x[1]}))
    variant_annotation$blacklist_region<- ifelse(pos %in% black_list_regions,"blacklist_region","")
    return(variant_annotation)
}


#' Annotate variants with mitochondrial homopolymer context
#'
#' @description
#' Flags each variant by the homopolymer region in which its position falls,
#' based on predefined mitochondrial homopolymer intervals.
#'
#' @details
#' The function loads the \code{mito_homopolymer} data frame, which must contain
#' columns \code{V2} (start), \code{V3} (end), and \code{V4} (region label).
#' For each variant position, it finds the first homopolymer interval that
#' contains that position and assigns its label; otherwise, \code{NA}.
#'
#' @param variant_annotation A data.frame or tibble containing at least a column
#'   named \code{Variants}, where each entry is of the form \code{"<pos>_<ref>_<alt>"}.
#'
#' @return The input \code{variant_annotation} with an added column:
#' \describe{
#'   \item{\code{homopolymer}}{Character; the label of the homopolymer region
#'     (from \code{mito_homopolymer\$V4}) containing the variant position, or
#'     \code{NA} if none.}
#' }
#'
#' @importFrom dplyr mutate
#' @export
annotate_variants_homopolymer <- function(variant_annotation){
    data(mito_homopolymer)   
    pos<- as.numeric(strsplit(variant_annotation$Variants, "_") %>% sapply(function(x){x[1]}))
    # for each pos, find the first interval that contains it
    assigned <- sapply(pos, function(pos) {
      idx <- which(mito_homopolymer$V2 <= pos & pos <= mito_homopolymer$V3)
      if (length(idx)) mito_homopolymer$V4[idx[1]] else NA
    })

    # assemble result
    pos_homopolymer <- data.frame(
      position = pos,
      V4        = assigned,
      stringsAsFactors = FALSE
    )
    variant_annotation <- variant_annotation %>% mutate(homopolymer = pos_homopolymer$V4)
    return(variant_annotation)
}


#' Parse mitochondrial variants into a tidy table
#'
#' @description
#' Internal helper to split variant strings of the form `"<pos>_<ref>_<mut>"`
#' into a data frame with columns for sample ID, chromosome, position, reference,
#' and alternate allele.
#'
#' @param variants Character vector of variant strings, e.g. `"123_A_G"`.
#' @param ID       Character scalar; sample identifier to assign to all rows.
#'
#' @return A data.frame with columns:
#'
#' @keywords internal
#' @noRd
make_mutation_table<-function(variants,ID=""){
    mito_mutations<-strsplit(variants,"_") %>% do.call(rbind,.) %>% as.data.frame()
    mito_mutations<-data.frame(sampleID=ID,chr="MT",mito_mutations)
    names(mito_mutations)[3:5]<-c("pos","ref","mut")
    mito_mutations$pos<-as.integer(mito_mutations$pos)
    return(mito_mutations)
}


#' Annotate variants with amino-acid change and predicted impact via dndscv
#'
#' @description
#' Runs dndscv on the mitochondrial variants to obtain gene-level
#' annotation of amino-acid changes and their predicted impact.
#'
#' @details
#' This internal function:
#' 1. Loads the precomputed mitochondrial CDS reference database from the
#'    redeemR package (`GRCH38_MT_bioMRT_output_refcds.rda`).
#' 2. Converts `Variants` ("<pos>_<ref>_<alt>") into the 5-column table
#'    required by dndscv using `make_mutation_table`.
#' 3. Runs `dndscv()` with permissive mutation limits to annotate all coding
#'    changes.
#' 4. Extracts and merges back the `gene`, `aachange`, and `impact` fields
#'    for each variant.
#'
#' @param variant_annotation A data.frame or tibble containing at least a column
#'   named `Variants` of the form "<pos>_<ref>_<alt>".
#'
#' @return The input `variant_annotation`, augmented with columns:
#' \describe{
#'   \item{gene}{Gene symbol where the variant falls.}
#'   \item{aachange}{Amino-acid change notation (e.g. p.W88C).}
#'   \item{impact}{Predicted functional impact category from dndscv.}
#' }
#'
#' @importFrom dndscv dndscv
#' @importFrom dplyr mutate select left_join any_of
#' @export
annotate_variants_aachange <- function(variant_annotation){
    require(dndscv)
    mitorefdb <- system.file("data", "GRCH38_MT_bioMRT_output_refcds.rda",package = "redeemR")
    mutations_all<- make_mutation_table(variant_annotation$Variants,"")
    mito_dndsout_all <- dndscv(mutations_all,refdb=mitorefdb,numcode = 2,max_muts_per_gene_per_sample = 5000,max_coding_muts_per_sample = 5000)
    mito_aa_annotation<- mito_dndsout_all$annotmuts %>% 
                         mutate(Variants= paste(pos, ref, mut, sep="_")) %>%
                         select(Variants, gene, aachange, impact)
    variant_annotation <- variant_annotation %>% select(-any_of(c("aachange", "impact"))) %>% left_join(mito_aa_annotation,by = "Variants")
}

#' Annotate variants with mitochondrial disease associations
#'
#' @description
#' Adds disease annotations from the mito_diseases dataset for each variant.
#'
#' @param variant_annotation A data.frame or tibble containing at least a column named \code{Variants}.
#'
#' @return The input \code{variant_annotation} with an added column:
#' \describe{\item{Disease}{Disease association(s) from \code{mito_diseases}, NA if none.}}
#' @importFrom dplyr select one_of left_join
#' @export
annotate_variants_mito_disease <- function(variant_annotation){
    data(mito_diseases)
    variant_annotation <- variant_annotation %>% 
                    select(-one_of(unique(c("Disease")))) %>%
                    left_join(mito_diseases[,c("ID","Disease")], by = c("Variants" = "ID"))
    return(variant_annotation)
}

#' Annotate hypermutable variants
#'
#' @description
#' Flags variants that fall into the hypermutable category based on cell-level 
#' mutation frequencies.
#'
#' @details
#' This function loads the \code{CellPCT.update} dataset and retrieves the set 
#' of hypermutable variants via \code{get_hype_v2(cut = 0.005)}. Those variant 
#' IDs are converted into the standard “<pos>_<ref>_<alt>” strings by 
#' \code{convert_variant()}, and any matches in the input table are labeled.
#'
#' @param variant_annotation A data.frame or tibble containing at least a column 
#'   named \code{Variants} whose entries are of the form “\<pos>_\<ref>_\<alt>”.
#'
#' @return The same \code{variant_annotation}, with an added column:
#' \describe{
#'   \item{\code{hyper_label}}{Character; \code{"hyper"} if the variant is  
#'     identified as hypermutable, otherwise \code{""}.}
#' }
#'
#' @importFrom dplyr mutate
#' @export
annotate_variants_hypermutable <- function(variant_annotation){
  data(CellPCT.update)
  hyper_v<-get_hype_v2(cut=0.005) %>% as.character %>% convert_variant()
  variant_annotation <- variant_annotation %>% mutate(hyper_label=ifelse(Variants %in% hyper_v, "hyper", ""))
  return(variant_annotation)
}


#' Apply full suite of variant annotations
#'
#' @description
#' Sequentially applies all annotation functions to a 
#' variant annotation data frame.
#'
#' @param variant_annotation A data.frame or tibble containing a column \\code{Variants}.
#'
#' @return The input \\code{variant_annotation}, with these annotations added in order:
#' \enumerate{
#'   \item Hypermutable labeling
#'   \item Population-frequency & haplogroup statistics
#'   \item Blacklist-region flagging
#'   \item Homopolymer-region annotation
#'   \item Amino-acid change & impact via dndscv
#'   \item Mitochondrial disease associations
#'   \item Base-change & transition/transversion classification
#' }
#'
#' @examples
#' df_annot <- annotate_all_variants(my_variant_df)
#' @importFrom dplyr %>%
#' @export
annotate_all_variants <- function(variant_annotation) {
  variant_annotation %>%
    annotate_variants_hypermutable() %>%
    annotate_variants_population_stats() %>%
    annotate_variants_blacklist() %>%
    annotate_variants_homopolymer() %>%
    annotate_variants_aachange() %>%
    annotate_variants_mito_disease() %>%
    Annotate_base_change()
}