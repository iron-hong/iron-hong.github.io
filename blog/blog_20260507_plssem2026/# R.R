# R
library(careless)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)

set.seed(123)

# --- 사용자 조정 파라미터 ---
add_synthetics <- TRUE        # 합성 불성실 응답 추가 여부
run_pct  <- 0.90              # max_run 임계값 퍼센타일 (상위)
sd_pct   <- 0.10              # sd 임계값 퍼센타일 (하위)
pt_pct   <- 0.10              # pt_corr 임계값 퍼센타일 (하위)
maha_pct <- 0.99              # mahalanobis 임계값 퍼센타일 (상위)
ridge_epsilon <- 1e-6         # 공분산 정규화 소량 값

# --- 데이터 준비 ---
df <- careless_dataset2
items <- df |> select(where(is.numeric))
item_names <- colnames(items)
n_items <- ncol(items)

# 선택적: 합성 불성실 응답(예시 3명)
if (add_synthetics) {
  straightliner <- as_tibble(setNames(rep(list(rep(4, n_items)), 1), NULL))
  colnames(straightliner) <- item_names
  random_responder <- as_tibble(setNames(list(sample(1:5, n_items, replace = TRUE)), NULL))
  colnames(random_responder) <- item_names
  longstring_responder <- {
    v <- c(rep(5, floor(n_items * 0.6)), sample(1:5, n_items - floor(n_items * 0.6), replace = TRUE))
    as_tibble(setNames(list(v), NULL))
  }
  colnames(longstring_responder) <- item_names
  items_aug <- bind_rows(items, straightliner, random_responder, longstring_responder)
} else {
  items_aug <- items
}

# --- 보조 함수 ---
max_run_length <- function(x) {
  r <- rle(x)
  if (length(r$lengths) == 0) return(0)
  max(r$lengths, na.rm = TRUE)
}
person_sd <- function(x) stats::sd(x, na.rm = TRUE)

# 공분산 정규화(특이행렬 방지)
center <- colMeans(items, na.rm = TRUE)
covmat <- cov(items, use = "pairwise.complete.obs")
covmat_reg <- covmat + diag(ridge_epsilon, ncol(covmat))

safe_maha <- function(x, center, covmat) {
  if (any(is.na(x))) return(NA_real_)
  tryCatch({
    mahalanobis(x, center, covmat)
  }, error = function(e) NA_real_)
}

# --- 지표 계산 ---
scores <- items_aug |>
  mutate(.row = row_number()) |>
  rowwise() |>
  mutate(
    vec = list(c_across(all_of(item_names))),
    max_run = max_run_length(vec),
    sd_resp = person_sd(vec),
    pt_corr = cor(vec, center, use = "pairwise.complete.obs"),
    maha = safe_maha(vec, center, covmat_reg)
  ) |>
  ungroup() |>
  select(.row, max_run, sd_resp, pt_corr, maha)

# --- 임계값 자동 결정(분포 기반) ---
thr_maxrun <- as.numeric(quantile(scores$max_run, run_pct, na.rm = TRUE))
thr_sd     <- as.numeric(quantile(scores$sd_resp, sd_pct, na.rm = TRUE))
thr_pt     <- as.numeric(quantile(scores$pt_corr, pt_pct, na.rm = TRUE))
thr_maha   <- as.numeric(quantile(scores$maha, maha_pct, na.rm = TRUE))

thresholds <- list(
  max_run = thr_maxrun,
  sd_low  = thr_sd,
  pt_corr = thr_pt,
  maha_high = thr_maha
)

# --- 플래그 규칙 적용 ---
flags <- scores |>
  mutate(
    flag_longstring = max_run >= thresholds$max_run,
    flag_invariant  = sd_resp <= thresholds$sd_low,
    flag_low_ptcorr = if_else(is.na(pt_corr), FALSE, pt_corr <= thresholds$pt_corr),
    flag_high_maha  = if_else(is.na(maha), FALSE, maha >= thresholds$maha_high),
    any_flag = flag_longstring | flag_invariant | flag_low_ptcorr | flag_high_maha
  )

# --- 시각화: 히스토그램(임계값선) 및 플래그별 상자그림 ---
p_maxrun_hist <- ggplot(flags, aes(x = max_run)) +
  geom_histogram(bins = 30, fill = "#2c7fb8", alpha = 0.6) +
  geom_vline(xintercept = thresholds$max_run, color = "red", linetype = "dashed") +
  labs(title = "Max run length distribution", x = "Max run", y = "Count")

p_sd_hist <- ggplot(flags, aes(x = sd_resp)) +
  geom_histogram(bins = 30, fill = "#7fc97f", alpha = 0.6) +
  geom_vline(xintercept = thresholds$sd_low, color = "red", linetype = "dashed") +
  labs(title = "Respondent SD distribution", x = "SD", y = "Count")

p_pt_hist <- ggplot(flags, aes(x = pt_corr)) +
  geom_histogram(bins = 30, fill = "#beaed4", alpha = 0.6) +
  geom_vline(xintercept = thresholds$pt_corr, color = "red", linetype = "dashed") +
  labs(title = "Person-total correlation distribution", x = "Pt-corr", y = "Count")

p_maha_hist <- ggplot(flags, aes(x = maha)) +
  geom_histogram(bins = 30, fill = "#fdc086", alpha = 0.6) +
  geom_vline(xintercept = thresholds$maha_high, color = "red", linetype = "dashed") +
  labs(title = "Mahalanobis distance distribution", x = "Mahalanobis", y = "Count")

# 플래그별 상자그림 (any_flag 기준)
flags_long <- flags |>
  pivot_longer(cols = c(max_run, sd_resp, pt_corr, maha),
               names_to = "metric", values_to = "value")

p_box_byflag <- ggplot(flags_long, aes(x = any_flag, y = value, fill = any_flag)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  facet_wrap(~ metric, scales = "free_y") +
  labs(x = "Flagged (any)", y = "Value") +
  theme(legend.position = "none")

# --- 결과 객체 반환 ---
result <- list(
  thresholds = thresholds,
  scores = scores,
  flags = flags,
  plots = list(
    maxrun_hist = p_maxrun_hist,
    sd_hist = p_sd_hist,
    pt_hist = p_pt_hist,
    maha_hist = p_maha_hist,
    box_by_flag = p_box_byflag
  )
)

result