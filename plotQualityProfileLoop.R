.plotQualityProfileLoop <- function(f_F = f_F, 
                                    f_R = f_R,
                                   sn = sn,
                                   plot_format = plot_format,
                                    fd = fd){
  fnFRs <- vector();
  for (i in seq(length(sn))){
    print(i);
    fnFRs[2*i -1] <- f_F[i];
    fnFRs[2*i] <- f_R[i];
  };
  n_plot_floor <- floor(length(sn) / 12);
  n_plot_remainder <- length(sn) %% 12;
  if (length(sn) <= 12){
    plotQualityProfile(fnFRs) -> p1; print(p1);
    ggsave(file.path(fd, 
                     paste0("sequence_quality.", plot_format)), 
           width = 20, height = 10);
  } else {
    for (i in seq(n_plot_floor)){
      plotQualityProfile(fnFRs[(i*12 - 11):(i*12)]) -> p2; print(p2);
      ggsave(file.path(fd,
                       paste0("sequence_quality_", i, ".", plot_format)),
             width = 20, height = 10)
    };
    if (n_plot_remainder > 0){
      plotQualityProfile(fnFRs[(n_plot_floor*12 + 1) : length(sn)]) -> p3; print(p3);
      ggsave(file.path(fd,
                       paste0("sequence_quality_", n_plot_floor+1, ".", plot_format)),
             width = 20, height = 10);
    };
  };
}
