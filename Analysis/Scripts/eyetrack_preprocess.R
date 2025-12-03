

pad_and_interpolate_blinks <- function(time, pupil, pad_ms = 100, sample_rate = 1000) {
  
  blink_mask = pupil == 0
  
  n_pad <- round(pad_ms / 1000 * sample_rate)
  
  padded_mask <- rep(FALSE, length(blink_mask))
  blink_indices <- which(blink_mask)
  
  if (length(blink_indices) == 0){
    return(pupil)
  }
  
  for (i in seq_along(blink_indices)) {
    idx <- blink_indices[i]
    start <- max(1, idx - n_pad)
    end   <- min(length(blink_mask), idx + n_pad)
    padded_mask[start:end] <- TRUE
  }
  
  interp_indices <- which(padded_mask)
  valid_indices <- which(!padded_mask)
  
  # Interpolate only if we have enough valid points
  if (length(valid_indices) > 1) {
    interp_pupil <- approx(x = time[valid_indices], y = pupil[valid_indices], xout = time)$y
    pupil[interp_indices] <- interp_pupil[interp_indices]
  }
  
  return(pupil)
}


interpolate_artifacts <- function(time, pupil, sample_rate = 1000,
                                  pad_ms = 100, artifact_sd_thresh = 5) {
  # Ensure numeric
  time <- as.numeric(time)
  pupil <- as.numeric(pupil)
  
  # Step 1: Detect abrupt changes (first derivative)
  pupil_diff <- c(0, diff(pupil))
  threshold <- artifact_sd_thresh * sd(pupil_diff, na.rm = TRUE)
  artifact_mask <- abs(pupil_diff) > threshold
  
  # Step 2: Pad around artifact samples
  n_pad <- round(pad_ms / 1000 * sample_rate)
  padded_mask <- rep(FALSE, length(artifact_mask))
  artifact_indices <- which(artifact_mask)
  
  for (idx in artifact_indices) {
    start <- max(1, idx - n_pad)
    end   <- min(length(artifact_mask), idx + n_pad)
    padded_mask[start:end] <- TRUE
  }
  
  # Step 3: Interpolate
  interp_indices <- which(padded_mask)
  valid_indices <- which(!padded_mask)
  
  interp_pupil <- pupil
  if (length(valid_indices) > 1 && !anyNA(time[valid_indices]) && !anyNA(pupil[valid_indices])) {
    interp_result <- approx(x = time[valid_indices], y = pupil[valid_indices], xout = time, rule = 2)$y
    interp_pupil[interp_indices] <- interp_result[interp_indices]
  } else {
    warning("Too few valid points for interpolation.")
  }
  return(interp_pupil)
  
  
}



