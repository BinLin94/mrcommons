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

## FAO's "stock"/"dairy"/"laying" (census) elements have year-specific reporting gaps -
## genuinely missing (NA) in the raw data, e.g. ITA poultry layers/broilers are NA in
## 2018, 2019, 2021, 2022 and 2024 (the latest available year) while the corresponding
## production stays positive throughout: a census gap, not a real population change.
## This is distinct from a real, reported zero (e.g. IRN pig stock, genuinely 0 every
## year since 1981, backed by production/slaughter also reported as 0 - a real,
## permanent stop). calcAnimalStocks() reads FAO_online with convert = "onlycorrect" and
## fills country gaps with NA (not 0) specifically so this function can tell the two
## apart from the raw data itself, rather than inferring it indirectly.
## For the categories derived by subtraction below ("other X = total stock - dairy/laying
## X"), an unfilled gap in either the total or the dairy/laying series sends that
## subtraction negative and the generic negative-value cleanup then wrongly zeroes out a
## real subcategory instead of just the corrupted residual - and even where the total is
## fine, a gap in just the dairy/laying series alone (e.g. ITA above) silently collapses
## MAgPIE's dairy/milk and egg-laying shares (calcLivestockDistribution's
## aggregateToMagpieKli) to zero even though the underlying species total is unaffected.
## For the categories used directly (swine, horses, ducks - no subtraction), the same gap
## just silently reads 0, with no negative value to flag it. All three cases are the same
## underlying data problem, so both the eight species stock totals and the five
## dairy/laying sub-population series that split them are gap-filled here before use,
## each paired against its own matching flow variable (meat production for the species
## totals, milk/egg production for the dairy/laying sub-populations):
##  - sandwiched (real values on both sides): linearly interpolated outright, since both
##    endpoints are hard evidence the population didn't hit zero in between.
##  - leading/trailing (run to the start/end of the series): filled if the paired meat
##    production/slaughter series is non-zero somewhere in the gap - evidence the
##    population is still there and it's the stock census specifically that's missing -
##    or if the gap is just a single year. Some countries report certain stock series
##    only periodically rather than annually (e.g. ITA poultry layers has recurring 1-2
##    year gaps in otherwise-stable years), so one missing year right at the current edge
##    of the data is far more likely to be "not published yet" than an instantaneous
##    change; longer gaps still require production evidence. Real zeros (like IRN pig
##    stock) never enter this function as a gap at all, since they are not NA.
##    Any gap still unresolved after this (no production evidence, longer than a single
##    year) is left as 0, e.g. IRN-style genuine, permanent stops would look the same if
##    they were ever reported as NA instead of a real 0.
toolFillStockGaps <- function(stock, production) {
  stockArr <- as.array(collapseNames(stock))[, , 1]
  prodArr  <- as.array(collapseNames(production))[, , 1]
  yrs      <- as.integer(gsub("^y", "", colnames(stockArr)))

  for (i in seq_len(nrow(stockArr))) {
    v <- stockArr[i, ]
    isGap <- is.na(v)
    if (!any(isGap) || all(isGap)) next # no gap, or never had any real data at all
    p <- prodArr[i, ]
    r <- rle(isGap)
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
        if (any(p[s:e] > 0, na.rm = TRUE) || (e - s + 1) == 1) v[s:e] <- v[s - 1]
      } else if (s == 1) {
        # leading
        if (any(p[s:e] > 0, na.rm = TRUE) || (e - s + 1) == 1) v[s:e] <- v[e + 1]
      }
    }
    stockArr[i, ] <- v
  }
  stockArr[is.na(stockArr)] <- 0 # any gap left unresolved (no evidence, longer than 1yr)

  out <- collapseNames(stock)
  out[, , ] <- as.numeric(stockArr)
  out
}

## Combine several raw FAO items (e.g. "Chickens" + "Other birds") into one series.
## Plain dimSums() without na.rm defaults to NA as soon as any one of the combined items
## is NA for a given country/year, even if another combined item has a real value there -
## which would wrongly flag every such country/year as a gap. na.rm = TRUE alone isn't
## enough either: base R (and dimSums) sums an all-NA set of inputs to 0, not NA, which
## would silently turn a genuine gap (e.g. a country that never reports "Other birds" at
## all, combined with a real gap year in "Chickens") into what looks like a real,
## reported zero - exactly the distinction this file exists to preserve. So: sum with
## na.rm = TRUE to let real values from other items carry a genuinely-missing one, then
## re-mark positions where every combined item was NA as NA again.
toolCombineItems <- function(x, items, dim) {
  sub   <- x[, , items]
  allNA <- Reduce(`&`, lapply(items, function(it) is.na(sub[, , it])))
  out   <- dimSums(sub, dim = dim, na.rm = TRUE)
  out[allNA] <- NA
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
  #
  # Read with convert = "onlycorrect" (skip convertFAO_online) and replicate its
  # remaining pieces manually - relative-value filtering, historical country mapping,
  # country-list completion - using fill = NA instead of the fill = 0 convertFAO_online
  # uses. That fill = 0 (shared by ~30 other calc functions reading FAO_online, so not
  # something to change there) conflates genuinely missing data with a real, reported
  # zero; toolFillStockGaps() above needs the distinction to only fill genuine gaps and
  # leave real zeros untouched. Any NA left over after gap-filling (categories not
  # passed through toolFillStockGaps, or countries never in FAO's data at all) is
  # cleaned up to 0 just before the return() at the end of this function.
  fao <- readSource("FAO_online", "LiveHead2024", convert = "onlycorrect")
  relativeDelete <- c("Yield_(100_g/ha)", "Yield_Carcass_Weight_(Hg/An)", "Yield_(100mg/An)",
                      "Yield_(No/An)", "Yield_(100_mg/An)", "Yield_(100_g/An)",
                      "Yield_Carcass_Weight_(100_g/An)", "Yield_Carcass_Weight_(0_1_g/An)",
                      "Yield_(100_g)")
  relativeDelete <- relativeDelete[relativeDelete %in% getItems(fao, dim = 3.2)]
  if (length(relativeDelete) > 0) fao <- fao[, , relativeDelete, invert = TRUE]

  # convert = "onlycorrect" skips convertFAO_online() entirely, which also skips its
  # handling of historical/composite FAO country codes (XET/XBL/XSD/XCN) that aren't in
  # madrat's default ISOhistorical mapping. Without this, toolCountryFill just drops them
  # as unknown codes - silently losing real data (e.g. Ethiopia's pre-1992 cattle stock is
  # entirely reported under "XET", not "ETH") instead of splitting/renaming it correctly.
  additionalMapping <- list()
  if (all(c("XET", "ETH", "ERI") %in% getItems(fao, dim = 1.1))) {
    additionalMapping <- append(additionalMapping, list(c("XET", "ETH", "y1992"), c("XET", "ERI", "y1992")))
  }
  if (all(c("XBL", "BEL", "LUX") %in% getItems(fao, dim = 1.1))) {
    additionalMapping <- append(additionalMapping, list(c("XBL", "BEL", "y1999"), c("XBL", "LUX", "y1999")))
  }
  if (all(c("XSD", "SSD", "SDN") %in% getItems(fao, dim = 1.1))) {
    additionalMapping <- append(additionalMapping, list(c("XSD", "SSD", "y2011"), c("XSD", "SDN", "y2011")))
  }
  if ("XCN" %in% getItems(fao, dim = 1.1)) {
    if ("CHN" %in% getItems(fao, dim = 1.1)) fao <- fao["CHN", , , invert = TRUE]
    getItems(fao, dim = 1)[getItems(fao, dim = 1) == "XCN"] <- "CHN"
  }

  fao <- toolISOhistorical(fao, overwrite = TRUE, additional_mapping = additionalMapping)
  fao <- toolCountryFill(fao, fill = NA, verbosity = 2)

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
  chickenStock <- toolFillStockGaps(toolCombineItems(liveHead, c("1057|Chickens", "1083|Other birds"), dim = 3.1),
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

  # dairy/laying sub-population series (used below to split each species' total stock
  # into dairy-or-layer vs other/broiler) have the same year-specific reporting gaps as
  # the total stocks above - gap-filled the same way, paired against milk/egg production
  # rather than meat production, since that is the flow variable that actually tracks
  # with these sub-populations
  dairyCowsStock  <- toolFillStockGaps(fao[, , "882|Raw milk of cattle.Milk_Animals_(An)"],
                                       fao[, , "882|Raw milk of cattle.Production_(t)"])
  dairyBufStock   <- toolFillStockGaps(fao[, , "951|Raw milk of buffalo.Milk_Animals_(An)"],
                                       fao[, , "951|Raw milk of buffalo.Production_(t)"])
  dairySheepStock <- toolFillStockGaps(fao[, , "982|Raw milk of sheep.Milk_Animals_(An)"],
                                       fao[, , "982|Raw milk of sheep.Production_(t)"])
  dairyGoatStock  <- toolFillStockGaps(fao[, , "1020|Raw milk of goats.Milk_Animals_(An)"],
                                       fao[, , "1020|Raw milk of goats.Production_(t)"])
  layerStock      <- toolFillStockGaps(
    toolCombineItems(fao, c("1062|Hen eggs in shell, fresh.Laying_(An)",
                            "1091|Eggs from other birds in shell, fresh, nec.Laying_(An)"), dim = 3),
    toolCombineItems(fao, c("1062|Hen eggs in shell, fresh.Production_(t)",
                            "1091|Eggs from other birds in shell, fresh, nec.Production_(t)"), dim = 3))

  # estimate numbers of animals for IPCC categories
  animals <- NULL

  # Dairy cows
  animals <- mbind(animals, setNames(
    collapseNames(dairyCowsStock), "dairy cows"))
  # Other cattle
  animals <- mbind(animals, setNames(
    cattleStock
    - setNames(animals[, , "dairy cows"], NULL),
    "other cattle"
  ))

  # Dairy Buffalo
  animals <- mbind(animals, setNames(
    collapseNames(dairyBufStock), "dairy buffalo"))
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
    collapseNames(dairySheepStock), "dairy sheep"))
  # Other sheep
  animals <- mbind(animals, setNames(
    sheepStock
    - setNames(animals[, , "dairy sheep"], NULL),
    "other sheep"
  ))

  # Dairy Goats
  animals <- mbind(animals, setNames(
    collapseNames(dairyGoatStock), "dairy goats"))
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
    toolCombineItems(liveHead, c("1126|Camels", "1157|Other camelids"), dim = 3.1)
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
    toolCombineItems(liveHead, c("1107|Asses", "1110|Mules and hinnies"), dim = 3.1),
    "mules and asses"
  ))

  # Poultry Layers
  animals <- mbind(animals, setNames(
    collapseNames(layerStock),
    "poultry layers"))
  # Broilers
  animals <- mbind(animals, setNames(
    chickenStock
    - setNames(animals[, , "poultry layers"], NULL),
    "broilers"
  ))

  # Turkey
  animals <- mbind(animals, setNames(
    toolCombineItems(liveHead, c("1072|Geese", "1079|Turkeys"), dim = 3.1),
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

  # categories not routed through toolFillStockGaps() (camels, turkey, mules/asses) keep
  # whatever NA survived toolCountryFill(fill = NA) - e.g. countries FAO never tracked
  # for that item at all. Clean up here rather than upstream, since upstream still needs
  # the NA/real-zero distinction.
  animals[is.na(animals)] <- 0

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
