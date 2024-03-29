% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nebula_re_rank.R
\name{nebula_re_rank}
\alias{nebula_re_rank}
\title{FUNCTION_TITLE}
\usage{
nebula_re_rank(
  nebula_name,
  top_n = 50,
  match_pattern = NULL,
  collate_factor = NA,
  only_gather_structure = F,
  reference_compound = NA,
  reference_ratio = 0.5,
  cluster_method = NA,
  csi_score_weight = 0.6,
  class_similarity_weight = 0.3,
  filter_via_classification = F,
  rt_set = NA,
  rt_weight = 0.1,
  rt_window = 1.5,
  revise_MCn_formula_set = F,
  revise_MCn_structure_set = F,
  ...
)
}
\arguments{
\item{nebula_name}{PARAM_DESCRIPTION}

\item{top_n}{PARAM_DESCRIPTION, Default: 10}

\item{match_pattern}{PARAM_DESCRIPTION, Default: c("precursorFormula")}

\item{collate_factor}{PARAM_DESCRIPTION, Default: 0.85}

\item{revise_MCn_formula_set}{PARAM_DESCRIPTION, Default: T}

\item{revise_MCn_structure_set}{PARAM_DESCRIPTION, Default: T}

\item{...}{PARAM_DESCRIPTION}
}
\value{
OUTPUT_DESCRIPTION
}
\description{
FUNCTION_DESCRIPTION
}
\details{
DETAILS
}
\examples{
\dontrun{
if(interactive()){
 #EXAMPLE1
 }
}
}
\seealso{
\code{\link[dplyr]{filter}}, \code{\link[dplyr]{reexports}}, \code{\link[dplyr]{arrange}}, \code{\link[dplyr]{mutate}}, \code{\link[dplyr]{select}}, \code{\link[dplyr]{distinct}}
\code{\link[pbapply]{pbapply}}
\code{\link[data.table]{rbindlist}}
\code{\link[tidyr]{separate}}
}
