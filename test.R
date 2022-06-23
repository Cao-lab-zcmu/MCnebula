## test
local({
  x <- 1
  local({
    x <<- 10
  })
})

