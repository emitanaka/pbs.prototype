test_that("pbs works", {
  set.seed(1)
  des <- design("Wheat yield & quality study") |>
    set_units(glasshouse = 4,
              block = nested_in(glasshouse, 4),
              pot = nested_in(block, 20)) |>
    set_trts(co2 = c("normal", "high"),
             temperature = c("normal", "high"),
             genotype = 20) |>
    allot_table(co2:temperature ~ glasshouse,
                genotype ~ pot)

  des |>
    set_rcrds(yield = pot, nutrition = pot) |>
    simulate_plant_experiment()

  des |>
    set_rcrds(yield = pot, nutrition = pot) |>
    simulate_plant_experiment(context = "the wheat used is durum wheat")


})
