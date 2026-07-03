#' Runtime Capabilities
#'
#' @description Reports runtime capabilities and platform information used by
#' \pkg{tsissm}. This is mainly useful for diagnosing platform-specific
#' numerical backends.
#'
#' @param force logical. If \code{TRUE}, re-run uncached runtime checks where
#' possible. User options such as \code{tsissm.rtmb_ad_eigen} still take
#' precedence.
#'
#' @return A named list with runtime capability information.
#' @export
#'
#' @examples
#' tsissm_capabilities()
tsissm_capabilities <- function(force = FALSE) {
    rtmb_ad_eigen_option <- getOption("tsissm.rtmb_ad_eigen", NULL)
    rtmb_ad_eigen <- rtmb_ad_eigen_available(force = force)
    rtmb_ad_eigen_source <- if (is.null(rtmb_ad_eigen_option)) {
        "runtime check"
    } else {
        "option"
    }
    list(
        rtmb_ad_eigen = rtmb_ad_eigen,
        rtmb_ad_eigen_source = rtmb_ad_eigen_source,
        rtmb_ad_eigen_option = rtmb_ad_eigen_option,
        RTMB = as.character(utils::packageVersion("RTMB")),
        TMB = as.character(utils::packageVersion("TMB")),
        R = paste(R.version$major, R.version$minor, sep = "."),
        platform = R.version$platform,
        os = Sys.info()[["sysname"]],
        release = Sys.info()[["release"]],
        machine = Sys.info()[["machine"]]
    )
}
