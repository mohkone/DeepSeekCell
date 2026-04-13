library(httr2)

api_key <- "sk-of-fShCqulKjvtVxdqcDnNfHBAjVVTEMNRtDBvGwOAlowmBWsxxRUqcrMikmVHHSCmw"

req <- request("https://api.ofox.ai/v1/chat/completions") |>
  req_headers(
    Authorization = paste("Bearer", api_key),
    "Content-Type" = "application/json"
  ) |>
  req_body_json(list(
    model = "openai/gpt-4o",
    messages = list(list(role = "user", content = "Say 'Hello'")),
    max_tokens = 50
  ))

resp <- req_perform(req)
data <- resp_body_json(resp)
cat(data$choices[[1]]$message$content)