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
##    change. Real zeros (like IRN pig stock) never enter this function as a gap at all,
##    since they are not NA.
##    Any gap still unresolved after this (no production evidence at all) is left as 0,
##    e.g. IRN-style genuine, permanent stops would look the same if they were ever
##    reported as NA instead of a real 0.
##
## A dairy/laying sub-population can also be genuinely 0 in a given year while its species
## total is a real positive value, with no within-series evidence to fill the
## sub-population itself (no sandwich, no matching milk/egg production in reach).
## toolImputeSubShare() below closes this case using that same country's own most recent
## real dairy/laying share (never another country's), clamped to [0, total] - see its own
## header for details.
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
        # sandwiched: real values on both sides, no evidence or distance limit needed
        v[s:e] <- approx(x = c(yrs[s - 1], yrs[e + 1]), y = c(v[s - 1], v[e + 1]), xout = yrs[s:e])$y
      } else if (e == length(v) || s == 1) {
        # trailing/leading: single-year gaps fill outright; longer gaps need production
        # evidence in each year
        trailing <- e == length(v)
        anchor   <- if (trailing) v[s - 1] else v[e + 1]
        supported <- (e - s + 1 == 1) | (!is.na(p[s:e]) & p[s:e] > 0)
        v[s:e][supported] <- anchor
      }
    }
    stockArr[i, ] <- v
  }
  stockArr[is.na(stockArr)] <- 0 # any gap left unresolved (no evidence) -> 0

  out <- collapseNames(stock)
  out[, , ] <- as.numeric(stockArr)
  out
}

## Impute a missing dairy/laying sub-population as (this country's most recent real
## dairy/laying share, last observation carried forward) x this year's total, clamped to
## [0, total] - never another country's share or a global average.
##
## `subRaw` (the sub-population as read from FAO, before toolFillStockGaps() ran) is
## required: the share anchor must be a real observation, not an estimate, and a real
## reported 0 (subRaw == 0) must never be overwritten - only a genuine gap (subRaw NA,
## still 0 in `sub`) is eligible.
toolImputeSubShare <- function(sub, total, subRaw) {
  subArr    <- as.array(collapseNames(sub))[, , 1]
  totalArr  <- as.array(collapseNames(total))[, , 1]
  subRawArr <- as.array(collapseNames(subRaw))[, , 1]

  for (i in seq_len(nrow(subArr))) {
    t   <- totalArr[i, ]
    raw <- subRawArr[i, ]
    # share anchor: only real, raw FAO observations count as evidence, never a value
    # toolFillStockGaps() already estimated (see header)
    validShare <- !is.na(raw) & raw > 0 & !is.na(t) & t > 0
    if (!any(validShare)) next # no real share evidence for this country/category at all
    shareSeries <- ifelse(validShare, raw / t, NA)
    # eligible only where the raw FAO reading was genuinely NA (not a real 0) and
    # toolFillStockGaps() still left it at 0 (no in-reach evidence to resolve it either)
    gap <- is.na(raw) & !is.na(t) & subArr[i, ] == 0 & t > 0
    if (!any(gap)) next

    for (yi in which(gap)) {
      priorValid <- which(validShare[seq_len(yi)])
      if (length(priorValid) == 0) next # no earlier real observation - leave unresolved
      share <- shareSeries[max(priorValid)]
      subArr[i, yi] <- min(max(share * t[yi], 0), t[yi])
    }
  }
  out <- collapseNames(sub)
  out[, , ] <- as.numeric(subArr)
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

  # toolISOhistorical() warns whenever the successor-country weight it looks up for some
  # unrelated item/element is NA in the split year - harmless (falls back to weight 0 for
  # just that item), but noisy given how many real NAs this file's fao object carries.
  # Silenced here; doesn't affect any values.
  fao <- suppressWarnings(toolISOhistorical(fao, overwrite = TRUE, additional_mapping = additionalMapping))
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

  # camels, turkey (geese+turkeys) and mules/asses each have a matching FAO meat-
  # production series, combined on both sides with toolCombineItems() before gap-filling,
  # same as chickenStock/layerStock below.
  camelStock <- toolFillStockGaps(
    toolCombineItems(liveHead, c("1126|Camels", "1157|Other camelids"), dim = 3.1),
    toolCombineItems(fao, c("1127|Meat of camels, fresh or chilled.Production_(t)",
                            "1158|Meat of other domestic camelids, fresh or chilled.Production_(t)"), dim = 3))
  turkeyStock <- toolFillStockGaps(
    toolCombineItems(liveHead, c("1072|Geese", "1079|Turkeys"), dim = 3.1),
    toolCombineItems(fao, c("1073|Meat of geese, fresh or chilled.Production_(t)",
                            "1080|Meat of turkeys, fresh or chilled.Production_(t)"), dim = 3))
  muleAssStock <- toolFillStockGaps(
    toolCombineItems(liveHead, c("1107|Asses", "1110|Mules and hinnies"), dim = 3.1),
    toolCombineItems(fao, c("1108|Meat of asses, fresh or chilled.Production_(t)",
                            "1111|Meat of mules, fresh or chilled.Production_(t)"), dim = 3))

  # dairy/laying sub-population series (used below to split each species' total stock
  # into dairy-or-layer vs other/broiler), gap-filled against milk/egg production rather
  # than meat production. Raw (pre-fill) copies are kept for toolImputeSubShare() below.
  dairyCowsRaw    <- fao[, , "882|Raw milk of cattle.Milk_Animals_(An)"]
  dairyCowsStock  <- toolFillStockGaps(dairyCowsRaw, fao[, , "882|Raw milk of cattle.Production_(t)"])
  dairyBufRaw     <- fao[, , "951|Raw milk of buffalo.Milk_Animals_(An)"]
  dairyBufStock   <- toolFillStockGaps(dairyBufRaw, fao[, , "951|Raw milk of buffalo.Production_(t)"])
  dairySheepRaw   <- fao[, , "982|Raw milk of sheep.Milk_Animals_(An)"]
  dairySheepStock <- toolFillStockGaps(dairySheepRaw, fao[, , "982|Raw milk of sheep.Production_(t)"])
  dairyGoatRaw    <- fao[, , "1020|Raw milk of goats.Milk_Animals_(An)"]
  dairyGoatStock  <- toolFillStockGaps(dairyGoatRaw, fao[, , "1020|Raw milk of goats.Production_(t)"])
  layerRaw        <- toolCombineItems(fao, c("1062|Hen eggs in shell, fresh.Laying_(An)",
                                             "1091|Eggs from other birds in shell, fresh, nec.Laying_(An)"), dim = 3)
  layerStock      <- toolFillStockGaps(
    layerRaw,
    toolCombineItems(fao, c("1062|Hen eggs in shell, fresh.Production_(t)",
                            "1091|Eggs from other birds in shell, fresh, nec.Production_(t)"), dim = 3))
  dairyCamelRaw   <- fao[, , "1130|Raw milk of camel.Milk_Animals_(An)"]
  dairyCamelStock <- toolFillStockGaps(dairyCamelRaw, fao[, , "1130|Raw milk of camel.Production_(t)"])

  # see toolImputeSubShare() header
  dairyCowsStock  <- toolImputeSubShare(dairyCowsStock, cattleStock, dairyCowsRaw)
  dairyBufStock   <- toolImputeSubShare(dairyBufStock, buffaloStock, dairyBufRaw)
  dairySheepStock <- toolImputeSubShare(dairySheepStock, sheepStock, dairySheepRaw)
  dairyGoatStock  <- toolImputeSubShare(dairyGoatStock, goatStock, dairyGoatRaw)
  layerStock      <- toolImputeSubShare(layerStock, chickenStock, layerRaw)
  dairyCamelStock <- toolImputeSubShare(dairyCamelStock, camelStock, dairyCamelRaw)

  # Externally-sourced corrections for specific country/category gaps that toolFillStockGaps()
  # leaves unresolved (no in-dataset production evidence) but a national statistics office
  # publishes a directly comparable (same species, same country, same unit) figure - found by
  # manually researching each of the ~39 country/category combinations still unresolved after
  # the fixes above. Deliberately narrow: only applied where a specific number was found and
  # cross-checked against an external source; only touches years still at the toolFillStockGaps
  # default of 0, never overwrites a year already resolved from FAO's own data. Countries/species
  # with only indirect ("the industry still exists, scale unconfirmed") evidence, or where the
  # external evidence was itself ambiguous (e.g. Palestine, where the gap coincides with the Gaza
  # war and FAO's PSE code combines Gaza with the still-reporting West Bank), are deliberately
  # left at the 0 default rather than guessed at. Full verification trail (including the entries
  # deliberately left unfilled) in the project's gap-audit record.
  applyExternalFill <- function(x, iso, years, value) {
    yrCols <- paste0("y", years)
    cur <- x[iso, yrCols, ]
    stillUnresolved <- as.numeric(cur) == 0
    cur[stillUnresolved] <- value
    x[iso, yrCols, ] <- cur
    x
  }

  # poultry layers (DPo)
  layerStock <- applyExternalFill(layerStock, "BGR", 2021:2024, 4392000) # Bulgaria Ministry of Agriculture and Food, 2024
  layerStock <- applyExternalFill(layerStock, "GLP", 2007:2024, 220155)  # DAAF Guadeloupe RA2020, 2020
  layerStock <- applyExternalFill(layerStock, "GUF", 2007:2024, 72000)   # DAAF Guyane Memento 2017, 2016
  layerStock <- applyExternalFill(layerStock, "IRL", 2021:2024, 3100000) # Eurostat egg-laying-hen statistics, 2022 (likely a slight
                                                                          # underestimate given IE's growing production trend since)
  layerStock <- applyExternalFill(layerStock, "MTQ", 2007:2024, 95000)   # last FAO value (2006) held constant - low confidence,
                                                                          # Agreste Martinique census notes an ongoing decline
  layerStock <- applyExternalFill(layerStock, "REU", 2007:2024, 508000)  # last FAO value (2006) held constant - medium confidence

  # ducks (Dk) - main species total, not a dairy/laying sub-population
  duckStock <- applyExternalFill(duckStock, "GUF", 2007:2024, 4700)   # DAAF Guyane Memento 2017, 2015 (700 canards a gaver +
                                                                        # 4,000 canards a rotir); last FAO value (30,000) was ~6x too high
  duckStock <- applyExternalFill(duckStock, "IRL", 2021:2024, 434000) # last FAO value (2020) held constant
  duckStock <- applyExternalFill(duckStock, "MTQ", 2002:2024, 590000) # last FAO value (2001) held constant - low confidence, decades old
  duckStock <- applyExternalFill(duckStock, "REU", 2007:2024, 520000) # last FAO value (2006) held constant

  # dairy cows (DCt)
  dairyCowsStock <- applyExternalFill(dairyCowsStock, "MTQ", 2007:2024, 96)    # DAAF Martinique Memento 2016, 2015;
                                                                                 # last FAO value (2,565, 2006) was ~27x too high
  dairyCowsStock <- applyExternalFill(dairyCowsStock, "REU", 2007:2024, 27280) # last FAO value (2006) held constant

  # dairy buffalo (DBf)
  dairyBufStock <- applyExternalFill(dairyBufStock, "ALB", 2020:2024, 29)    # last FAO value (2019) held constant; negligibly small
  dairyBufStock <- applyExternalFill(dairyBufStock, "PHL", 1991:2024, 30151) # Philippine Statistics Authority Dairy Situation Report,
                                                                               # 2023 - real growth since the 1990 baseline (8,000), so this
                                                                               # flat fill understates the 1990s-2000s and likely
                                                                               # overstates the 1991-1992 transition years slightly

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
    collapseNames(dairyCamelStock), "dairy camels"))
  # Other Camelids
  animals <- mbind(animals, setNames(
    camelStock
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
    muleAssStock,
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
    turkeyStock,
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

  # any NA remaining here is a country FAO never tracked for that item at all (dropped by
  # toolCountryFill(fill = NA) at read-in), not a within-series gap - toolFillStockGaps()
  # and toolImputeSubShare() already handled those. Cleaned up here rather than upstream,
  # since upstream still needs the NA/real-zero distinction.
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
