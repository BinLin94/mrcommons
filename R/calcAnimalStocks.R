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
#' @importFrom stats approx

## FAO's "stock" (census) elements have year-specific reporting gaps that read as a
## literal 0 even for large, active producers (e.g. DEU/ESP chicken stock is 0 in
## several years since 2018, while chicken meat production and slaughter counts stay
## normal throughout - a census gap, not a real population change). For the five
## categories derived by subtraction below ("other X = total stock - dairy/laying X"),
## an unfilled gap sends that subtraction negative and the generic negative-value
## cleanup then wrongly zeroes out a real subcategory instead of just the corrupted
## residual. For the three categories used directly (swine, horses, ducks - no
## subtraction), the same gap just silently reads 0, with no negative value to flag it.
## Both cases are the same underlying data problem, so all eight stock totals are
## gap-filled here before use:
##  - sandwiched (real values on both sides): linearly interpolated outright, since both
##    endpoints are hard evidence the population didn't hit zero in between.
##  - leading/trailing (run to the start/end of the series): only filled if the paired
##    meat production/slaughter series is non-zero somewhere in the gap - evidence the
##    population is still there and it's the stock census specifically that's missing.
##    Otherwise left as 0 (e.g. IRN pig stock 1981-2024, where production and slaughter
##    drop to 0 in the same year and never recover - a genuine, permanent stop, not a
##    reporting gap).
toolFillStockGaps <- function(stock, production) {
  stockArr <- as.array(collapseNames(stock))[, , 1]
  prodArr  <- as.array(collapseNames(production))[, , 1]
  yrs      <- as.integer(gsub("^y", "", colnames(stockArr)))

  for (i in seq_len(nrow(stockArr))) {
    v <- stockArr[i, ]
    isZero <- v == 0
    if (!any(isZero) || all(isZero)) next # no gap, or never had any stock at all
    p <- prodArr[i, ]
    r <- rle(isZero)
    ends   <- cumsum(r$lengths)
    starts <- ends - r$lengths + 1
    for (k in seq_along(r$lengths)) {
      if (!r$values[k]) next
      s <- starts[k]; e <- ends[k]
      if (s > 1 && e < length(v)) {
        # sandwiched: real values on both sides
        v[s:e] <- approx(x = c(yrs[s - 1], yrs[e + 1]), y = c(v[s - 1], v[e + 1]), xout = yrs[s:e])$y
      } else if (e == length(v)) {
        # trailing
        if (any(p[s:e] > 0, na.rm = TRUE)) v[s:e] <- v[s - 1]
      } else if (s == 1) {
        # leading
        if (any(p[s:e] > 0, na.rm = TRUE)) v[s:e] <- v[e + 1]
      }
    }
    stockArr[i, ] <- v
  }

  out <- collapseNames(stock)
  out[, , ] <- as.numeric(stockArr)
  out
}

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

  # gap-filled versions of the raw stock totals feeding the subtractions below -
  # see toolFillStockGaps() above for why this is needed
  cattleStock  <- toolFillStockGaps(liveHead[, , "866|Cattle"],
                                    fao[, , "867|Meat of cattle with the bone, fresh or chilled.Production_(t)"])
  buffaloStock <- toolFillStockGaps(liveHead[, , "946|Buffalo"],
                                    fao[, , "947|Meat of buffalo, fresh or chilled.Production_(t)"])
  sheepStock   <- toolFillStockGaps(liveHead[, , "976|Sheep"],
                                    fao[, , "977|Meat of sheep, fresh or chilled.Production_(t)"])
  goatStock    <- toolFillStockGaps(liveHead[, , "1016|Goats"],
                                    fao[, , "1017|Meat of goat, fresh or chilled.Production_(t)"])
  chickenStock <- toolFillStockGaps(dimSums(liveHead[, , c("1057|Chickens", "1083|Other birds")], dim = 3.1),
                                    fao[, , "1058|Meat of chickens, fresh or chilled.Production_(t)"])
  # these three are used directly (no subtraction), so a stock gap never goes negative
  # and the cleanup at the end of this function never flags it - it just silently reads
  # 0. Gap-filled here too so the underlying data quality issue is fixed the same way.
  swineStock   <- toolFillStockGaps(liveHead[, , "1034|Swine / pigs"],
                                    fao[, , "1035|Meat of pig with the bone, fresh or chilled.Production_(t)"])
  horseStock   <- toolFillStockGaps(liveHead[, , "1096|Horses"],
                                    fao[, , "1097|Horse meat, fresh or chilled.Production_(t)"])
  duckStock    <- toolFillStockGaps(liveHead[, , "1068|Ducks"],
                                    fao[, , "1069|Meat of ducks, fresh or chilled.Production_(t)"])

  # estimate numbers of animals for IPCC categories
  animals <- NULL

  # Dairy cows
  animals <- mbind(animals, setNames(
    collapseNames(fao[, , "882|Raw milk of cattle.Milk_Animals_(An)"]), "dairy cows"))
  # Other cattle
  animals <- mbind(animals, setNames(
    cattleStock
    - setNames(animals[, , "dairy cows"], NULL),
    "other cattle"
  ))

  # Dairy Buffalo
  animals <- mbind(animals, setNames(
    collapseNames(fao[, , "951|Raw milk of buffalo.Milk_Animals_(An)"]), "dairy buffalo"))
  # Other buffalo
  animals <- mbind(animals, setNames(
    buffaloStock
    - setNames(animals[, , "dairy buffalo"], NULL), "other buffalo"))

  # Market Swine
  animals <- mbind(animals, setNames(
    swineStock * marketSwineShare,
    "market swine"
  ))

  # Breeding Swine
  animals <- mbind(animals, setNames(
    swineStock * (1 - marketSwineShare),
    "breeding swine"
  ))

  # Dairy Sheep
  animals <- mbind(animals, setNames(
    collapseNames(fao[, , "982|Raw milk of sheep.Milk_Animals_(An)"]), "dairy sheep"))
  # Other sheep
  animals <- mbind(animals, setNames(
    sheepStock
    - setNames(animals[, , "dairy sheep"], NULL),
    "other sheep"
  ))

  # Dairy Goats
  animals <- mbind(animals, setNames(
    collapseNames(fao[, , "1020|Raw milk of goats.Milk_Animals_(An)"]), "dairy goats"))
  # Other goats
  animals <- mbind(animals, setNames(
    goatStock
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
    horseStock,
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
    chickenStock
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
    duckStock,
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
