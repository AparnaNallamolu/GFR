#' Cross-validation with Random partitions
#'
#' @name CV.RandomPart
#' @rdname CV.RandomPart
#' @keywords internal
#' @export
#' @importFrom BMTME CV.RandomPart
#' @usage CV.RandomPart(data, Npartitions = 10,  PTesting = 0.2, set_seed = 123)
NULL


#' Cross-validation with K Folds
#'
#' @name CV.KFold
#' @rdname CV.KFold
#' @keywords internal
#' @export
#' @importFrom BMTME CV.KFold
#' @usage CV.KFold(data, DataSetID = 'Line', K = 5, set_seed =  123)
NULL

#' Cross-validation with K Folds
#'
#' @name CV.KFold
#' @rdname CV.KFold
#' @keywords internal
#' @export
#' @importFrom BMTME CV.KFold
#' @usage CV.KFold(data, DataSetID = 'Line', K = 5, set_seed =  123)
NULL

#' Cross-validation with stratified samples
#'
#' @name CV.StratifiedByFrac
#' @rdname CV.StratifiedByFrac
#' @keywords internal
#' @export
#' @importFrom BMTME CV.StratifiedByFrac
#' @usage CV.StratifiedByFrac(DataSet, NSamples = 10, fracTesting = 0.1, replace = FALSE, set_seed = NULL)
NULL

#' Cross-Validation with stratified samples
#'
#' @name CV.Stratified
#' @rdname CV.Stratified
#' @keywords internal
#' @export
#' @importFrom BMTME CV.Stratified
#' @usage CV.Stratified(DataSet, NSamples = 10, nTesting = 10, replace = FALSE, set_seed = NULL)
NULL
