#' Describe context
#'
#' `describe_context()` describes the context for the given edibble experiment
#' in words. The `describe_scenario()` takes the former context and uses a LLM
#' to create an embellished version of the context.
#'
#' @param .edibble An edibble design.
#'
#' @export
describe_context <- function(.edibble) {
  not_edibble(.edibble)
  prov <- activate_provenance(.edibble)
  prov$get_title()

  desc <- paste0("Title: ", prov$get_title())
  if(length(units <- prov$fct_names(role = "edbl_unit"))) {
    desc <- c(desc, paste0("Units: ", paste0(units, collapse = ", ")))
  }
  if(length(trts <- prov$fct_names(role = "edbl_trt"))) {
    desc <- c(desc, paste0("Treatments: ", paste0(trts, collapse = ", ")))
  }
  if(length(rcrds <- prov$fct_names(role = "edbl_rcrd"))) {
    desc <- c(desc, paste0("Records: ", paste0(rcrds, collapse = ", ")))
  }
  if(nrow(edges <- prov$fct_edges)) {
    verbs <- map_chr(edges$type, function(type) switch(type,
                                                       nest = "nested in",
                                                       allot = "assigned to",
                                                       record = "measured on",
                                                       "related to"))

    desc <- c(desc, paste0("Mappings: ", paste0(paste(edges$var_from,
                                                      "is",
                                                      verbs,
                                                      edges$var_to),
                                                collapse = ", ")))
  }

  define_context(paste0(desc, collapse = "\n"),
                 design = .edibble)
}

#' @rdname describe_context
#' @export
describe_scenario <- function(.edibble, chat = ellmer::chat_openai()) {
  context <- describe_context(.edibble)
  chat$set_system_prompt("Given the structure of the experiment, describe the experiment.")
  define_context(chat$chat(context),
                 design = .edibble,
                 chat = chat)
}

#' Check if it is an edibble context
#'
#' An edibble context would be a text with the edibble design
#' in the attribute and optionally a chat object.
#'
#' @param x An object to test.
#'
#' @export
is_edibble_context <- function(x) {
  inherits(x, "edbl_context")
}


define_context <- function(x, design, chat = NULL) {
  x <- glue::glue(x)
  structure(x,
            class = c("edbl_context", class(x)),
            design = edbl_design(design),
            chat = chat)
}

extract_chat <- function(x) {
  attr(x, "chat")
}
