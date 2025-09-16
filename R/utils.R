

`%||%` <- function(x, y) {
  if(is.null(x)) y else x
}

not_edibble <- function(x) {
  if (!is_edibble(x)) {
    rlang::abort(sprintf("%s is not an edibble.", deparse(substitute(x))))
  }
}

glue_c <- function(content, last = " and ", .content_envir = rlang::caller_env()) {
  glue::glue(content, .transformer = collapse_transformer("[*]$", last = last), .envir = .content_envir)
}

make_robj_name <- function(text) {
  make.names(text)
}

# as inspired by the glue doc
collapse_transformer <- function(regex, last) {
  function(text, envir) {
    collapse <- grepl(regex, text)
    if (collapse) {
      text <- sub(regex, "", text)
    }
    res <- glue::identity_transformer(text, envir)
    if (collapse) {
      glue::glue_collapse(res, sep = ", ", last = last)
    } else {
      res
    }
  }
}


