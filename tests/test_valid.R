#### ValidCharacter ####

test_that("validCharacter", {
  test.value <- as.character(1)
  validCharacter(test.value, validLength = 1)
  test.value <- as.character(1:2)
  validCharacter(test.value, validLength = 2)
  test.value <- as.character(1:2)
  validCharacter(test.value, validLength = NULL)
  
  test.value <- 1:2
  expect_error(validCharacter(test.value, validLength = 1))
  expect_error(validCharacter(test.value, validLength = 2))

  test.value <- c(T,F)
  validCharacter(test.value, validLength = 2, validValues = "character_or_logical")
  expect_error(validCharacter(test.value, validLength = 2, validValues = "character"))
  test.value <- c("T","F")
  validCharacter(test.value, validLength = 2, validValues = "character")
  
  test.value <- NULL
  validCharacter(test.value, validLength = NULL, refuse.NULL = FALSE)
  expect_error(validCharacter(test.value, validLength = NULL, refuse.NULL = TRUE))
  
  test.value <- rep("1",2)
  validCharacter(test.value, validLength = NULL, refuse.duplicates = FALSE)
  expect_error(validCharacter(test.value, validLength = NULL, refuse.duplicates = TRUE))
})

#### validClass ####

#### validDimension ####

#### validInteger ####

#### validNames ####

#### validNumeric ####

#### validPath ####

