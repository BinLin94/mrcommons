#' @title calcAnimalStocks
#' @description calculates stocks of animals of different categories.
#'
#' @param grouping IPCC: Animal grouping of IPCC Guidelines
#'
#' @return List of magpie objects with results on country level, weight on country level, unit and description.
#' @author Benjamin Leon Bodirsky, Bin Lin
#' @seealso
#' [calcExcretionIPCC()],
#' [readIPCC()]
#' @examples
#' \dontrun{
#' calcOutput("AnimalStocks")
#' }
#'
calcAnimalStocks <- function(grouping = "IPCC") {
  if (grouping != "IPCC") {
    stop("so far only IPCC categories implemented.")
  }

  marketSwineShare <- 0.9 # table 10.19

  # FAO merged LiveHead/LivePrim into Production_Crops_Livestock in 2024.
  # Use LiveHead2024 which reads the new merged file; item names changed to
  # e.g. "882|Raw milk of cattle" and elements to "Milk_Animals_(An)".
  fao      <- readSource("FAO_online", "LiveHead2024")
  liveHead <- dimSums(fao, dim = "ElementShort")

  # estimate numbers of animals for IPCC categories
  animals <- NULL

  # Dairy cows
  animals <- mbind(animals, setNames(
    collapseNames(fao[, , "882|Raw milk of cattle.Milk_Animals_(An)"]), "dairy cows"))
  # Other cattle
  animals <- mbind(animals, setNames(
    collapseNames(liveHead[, , "866|Cattle"])
    - setNames(animals[, , "dairy cows"], NULL),
    "other cattle"
  ))

  # Dairy Buffalo
  animals <- mbind(animals, setNames(
    collapseNames(fao[, , "951|Raw milk of buffalo.Milk_Animals_(An)"]), "dairy buffalo"))
  # Other buffalo
  animals <- mbind(animals, setNames(
    collapseNames(liveHead[, , "946|Buffalo"])
    - setNames(animals[, , "dairy buffalo"], NULL), "other buffalo"))

  # Market Swine
  animals <- mbind(animals, setNames(
    collapseNames(liveHead[, , "1034|Swine / pigs"]) * marketSwineShare,
    "market swine"
  ))

  # Breeding Swine
  animals <- mbind(animals, setNames(
    collapseNames(liveHead[, , "1034|Swine / pigs"]) * (1 - marketSwineShare),
    "breeding swine"
  ))

  # Dairy Sheep
  animals <- mbind(animals, setNames(
    collapseNames(fao[, , "982|Raw milk of sheep.Milk_Animals_(An)"]), "dairy sheep"))
  # Other sheep
  animals <- mbind(animals, setNames(
    collapseNames(liveHead[, , "976|Sheep"])
    - setNames(animals[, , "dairy sheep"], NULL),
    "other sheep"
  ))

  # Dairy Goats
  animals <- mbind(animals, setNames(
    collapseNames(fao[, , "1020|Raw milk of goats.Milk_Animals_(An)"]), "dairy goats"))
  # Other goats
  animals <- mbind(animals, setNames(
    collapseNames(liveHead[, , "1016|Goats"])
    - setNames(animals[, , "dairy goats"], NULL),
    "other goats"
  ))

  # Dairy Camels
  animals <- mbind(animals, setNames(
    collapseNames(fao[, , "1130|Raw milk of camel.Milk_Animals_(An)"]), "dairy camels"))
  # Other Camelids
  animals <- mbind(animals, setNames(
    dimSums(liveHead[, , c("1126|Camels", "1157|Other camelids")], dim = 3.1)
    - setNames(animals[, , "dairy camels"], NULL),
    "other camels"
  ))

  # Horses
  animals <- mbind(animals, setNames(
    collapseNames(liveHead[, , "1096|Horses"]),
    "horses"
  ))

  # Mules and Asses
  animals <- mbind(animals, setNames(
    dimSums(liveHead[, , c("1107|Asses", "1110|Mules and hinnies")], dim = 3.1),
    "mules and asses"
  ))

  # Poultry Layers
  animals <- mbind(animals, setNames(
    dimSums(fao[, , c("1062|Hen eggs in shell, fresh.Laying_(An)",
                      "1091|Eggs from other birds in shell, fresh, nec.Laying_(An)")], dim = 3),
    "poultry layers"))
  # Broilers
  animals <- mbind(animals, setNames(
    dimSums(liveHead[, , c("1057|Chickens", "1083|Other birds")], dim = 3.1)
    - setNames(animals[, , "poultry layers"], NULL),
    "broilers"
  ))

  # Turkey
  animals <- mbind(animals, setNames(
    dimSums(liveHead[, , c("1072|Geese", "1079|Turkeys")], dim = 3.1),
    "turkey"
  ))

  # Ducks
  animals <- mbind(animals, setNames(
    collapseNames(liveHead[, , "1068|Ducks"]),
    "ducks"
  ))

  # ignore
  # "1140|Rabbits and hares","1150|Rodents, other"
  animals <- animals / 1000000

  # sort according to n_rate animal categories


  # remove all negative values
  remove <- which(animals < 0)

  if (length(remove) > 0) {
    vcat(2, paste0(length(remove), " negative values removed"))
    animals[remove] <- 0
  }

  return(list(
    x = animals,
    weight = NULL,
    unit = "Million animals",
    description = "Animal stocks, for laying hens and dairy cattle producing animals",
    min = 0
  ))
}
