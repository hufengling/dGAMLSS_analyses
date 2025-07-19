generate_penalty_matrix <- function(local_object,
                                      smooth_penalty_list,
                                      smooth_index_list) {
  penalty_matrix_list <- list()
  if ("mu" %in% local_object$parameters) {
    if (("mu" %in% names(smooth_penalty_list) & !("mu" %in% names(smooth_index_list))) |
        (!("mu" %in% names(smooth_penalty_list)) & "mu" %in% names(smooth_index_list))) {
      stop("Both smooth_penalty_list and smooth_index_list must be named lists with the same penalized parameters as the names.")
    }

    if ("mu" %in% names(smooth_penalty_list) & "mu" %in% names(smooth_index_list)) {
      tmp_matrix <- smooth_penalty_list$mu
      tmp_index <- smooth_index_list$mu
      if ((sum(tmp_index) != nrow(tmp_matrix)) | !all(tmp_index == 0 | tmp_index == 1)) {
        stop("smooth_index_list must be a list of vectors of 0s and 1s indicating which coefficient indices correspond to the spline. Sum of smooth_index_list vectors must be equal to sides of corresponding smooth_penalty_list matrices.")
      }
      tmp_combined <- tmp_index %*% t(tmp_index)
      tmp_combined[tmp_combined == 1] <- tmp_matrix
      penalty_matrix_list$mu <- tmp_combined
    }
  }

  if ("sigma" %in% local_object$parameters) {
    if (("sigma" %in% names(smooth_penalty_list) & !("sigma" %in% names(smooth_index_list))) |
        (!("sigma" %in% names(smooth_penalty_list)) & ("sigma" %in% names(smooth_index_list)))) {
      stop("Both smooth_penalty_list and smooth_index_list must be named lists with the same penalized parameters as the names.")
    }

    if ("sigma" %in% names(smooth_penalty_list) & "sigma" %in% names(smooth_index_list)) {
      tmp_matrix <- smooth_penalty_list$sigma
      tmp_index <- smooth_index_list$sigma
      if ((sum(tmp_index) != nrow(tmp_matrix)) | !all(tmp_index == 0 | tmp_index == 1)) {
        stop("smooth_index_list must be a list of vectors of 0s and 1s indicating which coefficient indices correspond to the spline. Sum of smooth_index_list vectors must be equal to sides of corresponding smooth_penalty_list matrices.")
      }
      tmp_combined <- tmp_index %*% t(tmp_index)
      tmp_combined[tmp_combined == 1] <- tmp_matrix
      penalty_matrix_list$sigma <- tmp_combined
    }
  }

  if ("nu" %in% local_object$parameters) {
    if (("nu" %in% names(smooth_penalty_list) & !("nu" %in% names(smooth_index_list))) |
        (!("nu" %in% names(smooth_penalty_list)) & "nu" %in% names(smooth_index_list))) {
      stop("Both smooth_penalty_list and smooth_index_list must be named lists with the same penalized parameters as the names.")
    }

    if ("nu" %in% names(smooth_penalty_list) & "nu" %in% names(smooth_index_list)) {
      tmp_matrix <- smooth_penalty_list$nu
      tmp_index <- smooth_index_list$nu
      if ((sum(tmp_index) != nrow(tmp_matrix)) | !all(tmp_index == 0 | tmp_index == 1)) {
        stop("smooth_index_list must be a list of vectors of 0s and 1s indicating which coefficient indices correspond to the spline. Sum of smooth_index_list vectors must be equal to sides of corresponding smooth_penalty_list matrices.")
      }
      tmp_combined <- tmp_index %*% t(tmp_index)
      tmp_combined[tmp_combined == 1] <- tmp_matrix
      penalty_matrix_list$nu <- tmp_combined
    }
  }

  if ("tau" %in% local_object$parameters) {
    if (("tau" %in% names(smooth_penalty_list) & !("tau" %in% names(smooth_index_list))) |
        (!("tau" %in% names(smooth_penalty_list)) & "tau" %in% names(smooth_index_list))) {
      stop("Both smooth_penalty_list and smooth_index_list must be named lists with the same penalized parameters as the names.")
    }

    if ("tau" %in% names(smooth_penalty_list) & "tau" %in% names(smooth_index_list)) {
      tmp_matrix <- smooth_penalty_list$tau
      tmp_index <- smooth_index_list$tau
      if ((sum(tmp_index) != nrow(tmp_matrix)) | !all(tmp_index == 0 | tmp_index == 1)) {
        stop("smooth_index_list must be a list of vectors of 0s and 1s indicating which coefficient indices correspond to the spline. Sum of smooth_index_list vectors must be equal to sides of corresponding smooth_penalty_list matrices.")
      }
      tmp_combined <- tmp_index %*% t(tmp_index)
      tmp_combined[tmp_combined == 1] <- tmp_matrix
      penalty_matrix_list$tau <- tmp_combined
    }
  }

  penalty_matrix_list$smooth_index_list <- smooth_index_list
  penalty_matrix_list
}
