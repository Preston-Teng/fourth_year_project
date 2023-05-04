
#———load and download all the relevant packages———

library(tidyverse)
library(litsearchr)
library(igraph)
library(here)

#———importing naive search results———
#———naive search done on WOS + Scopus + BASE———

naive_import <-
  import_results(
    directory = here("naive_results"),
    verbose = TRUE)

#———remove duplicates———

naive_results <- remove_duplicates(naive_import, field = "title", method = "exact")

#———establish stop words———

list_of_stopwords <- c("phytoplankton",
                       "coastal",
                       "coast",
                       "lake",
                       "continental",
                       "fish",
                       "fishes",
                       "estuary",
                       "straits",
                       "ice",
                       "bacteria",
                       "isolated",
                       "krill",
                       "euphausiids",
                       "bank"
)

all_stopwords <- c(get_stopwords("English"), list_of_stopwords)

#———extract keywords———

my_text <- paste(naive_results$title,
                 naive_results$abstract,
                 naive_results$keywords)

raked_keywords <- extract_terms(
  text = my_text,
  method = "fakerake",
  min_freq = 2,
  ngrams = TRUE,
  min_n = 2,
  language = "English",
  stopwords = all_stopwords
)


#———building a co-occurrence network———

naive_dfm <- create_dfm(
  elements = my_text,
  features = raked_keywords
)

naive_graph <- create_network(
  search_dfm = as.matrix(naive_dfm),
  min_studies = 2,
  min_occ = 2
)

#plot_network(graph = naive_graph,
#            graphcolor = "#006e6d",
#           color_gradient = TRUE)
#---> not sure how to plot the net work
#---> should check

#———identifying keyword importance/ strength———

strengths <- sort(strength(naive_graph), decreasing=TRUE)

# to see what the top 10 terms are head(strengths, 10)

# if we want to plot the occurence terms --> plot(strengths, type="l", las=1)

plot(strengths, type="l", las=1)

term_strengths <- data.frame(term = names(strengths),
                             strength = strengths,
                             row.names = NULL) %>%
  mutate(rank = rank(strength, ties.method = "min")) %>%
  arrange(strength)

cutoff <- find_cutoff(
  naive_graph,
  method = "cumulative",
  percent = .30, #what percentage we want to use --> in the example it is 0.5
  imp_method = "strength"
)

#———visualising the co-occurrence network———

ggraph(naive_graph, layout="stress") +
  coord_fixed() +
  expand_limits(x=c(-3, 3)) +
  geom_edge_link(aes(alpha=weight)) +
  geom_node_point(shape="circle filled", fill="white") +
  geom_node_text(aes(label=name), hjust="outward", check_overlap=TRUE) +
  guides(edge_alpha=FALSE)

ggplot(term_strengths, aes(x=rank, y=strength, label=term)) +
  geom_line() +
  geom_point() +
  geom_text(data = filter(term_strengths, rank > 5),
            hjust = "right",
            nudge_y = 20,
            check_overlap = TRUE) +
  geom_hline(yintercept = cutoff, linetype = "dashed")

#———creating a reduced network———

reduced_graph <- reduce_graph(naive_graph, cutoff_strength = cutoff)

search_terms <- get_keywords(reduced_graph)

head(search_terms, 30)

#———group terms———

write.csv(search_terms, "./suggested_search_terms.csv")

grouped_terms <- read.csv("./suggested_search_terms_final.csv")

head(grouped_terms)

population <- grouped_terms$term[grep(pattern="population", grouped_terms$group)]

method <- grouped_terms$term[grep(pattern="method", grouped_terms$group)]

outcome <- grouped_terms$term[grep(pattern="outcome", grouped_terms$group)]

#———appending additional terms———

population <- append(population, c("pacific ocean",
                                   "arctic ocean",
                                   "atlantic ocean",
                                   "indian ocean",
                                   "north atlantic",
                                   "south atlantic",
                                   "antarctic ocean",
                                   "mesozooplankton",
                                   "zooplankton community"))

exposure <- append(method, c("cruise",
                             "station"))

outcome <- append(outcome, c("biomass integrated",
                             "biomass",
                             "zooplankton biomass",
                             "abundance",
                             "abundance integrated",
))

#———writing a boolean search———

full_search <- write_search(
  list(population, method, outcome),
  closure = "none",
  languages = "English",
  stemming = TRUE,
  exactphrase = TRUE,
  writesearch = FALSE,
  verbose = TRUE
)

print(full_search)

write.csv(full_search, “./final_full_search_terms.csv”)


