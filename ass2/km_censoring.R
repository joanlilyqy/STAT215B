N = 100
end.of.study = 5


gen.km = function(event.times, was.death) {
  
  event.ordering = order(event.times)
  
  s.hat.y = c(1)
  
  for (i in 1:N) {
    e.id = event.ordering[i]
    if (event.times[e.id] > end.of.study)
      break
    if (was.death[e.id]) {
      ni = N - i + 1  # number at risk
      next.factor = (ni - 1) / ni
      prev.factor = s.hat.y[length(s.hat.y)]
      s.hat.y = c(s.hat.y, prev.factor * next.factor)
    }
  }
  
  observed.death.times = sort(event.times[was.death])
  s.hat.x = observed.death.times[observed.death.times < end.of.study]

  # 1. return a function (could be your own, but it should be 
  #     fast enough to finish in seconds)
  # 2. the function return should work for any t
  # 3. no need to return anything else
  stepfun(s.hat.x, s.hat.y)
}


plot.km = function() {
  plot(km, 
       xlim=c(0, 5),
       ylim=c(0,1), 
       do.points=F,
       xlab="study duration (years)", 
       ylab="survival rate", 
       main="", 
       col="red")
  
  # this is the true curve, not the some uncensored KM estimator
  curve(exp(-(x/2)^2), 
        from=0, 
        to=5, 
        add=T, 
        col = "black")
  
  legend("topright", 
         lty=1, 
         col=c("black", "red"), 
         legend=c("true", "Kaplan-Meier"),
         cex=.75)
}



death.times = rweibull(N, shape=2, scale=2)


# hw 2, question 3: no censoring

was.death = rep(T, N)
km = gen.km(death.times, was.death)
plot.km()


# hw 2, question 4: independent censoring

censor.times = rexp(N, rate=.1)
was.death = death.times < censor.times
#event times are the x coordinates, not the sorted death times
event.times = sapply(1:N, function(i) min(death.times[i], censor.times[i]))
km = gen.km(event.times, was.death)
plot.km()



# hw 2, question 5: dependent censoring

# subjects who lived less than 2 years were
# less likely to quit the study.

# extreme censoring through thresholding may shift the KM curve up/out...

get.censored = function(i)
  if (death.times[i] >= 2) 2 else censor.times[i]
alt.censor.times = sapply(1:N, get.censored)
row.min = function(i) min(death.times[i], alt.censor.times[i])
alt.event.times = sapply(1:N, row.min)
alt.was.death = death.times < alt.censor.times

alt.km = gen.km(alt.event.times, alt.was.death)
plot.km()
lines(alt.km, col="blue", do.points=F, xlim=c(0,5))
#Why do the lines deviate between 0 and 2?


# but the question is about moderate censoring though exponential filtering.
# in this senario, subjects appear to be dying sooner.

# the min of exponentials is itself exponential with rate \lambda_1 + \lambda_2
# reusing the previous exp reduces noise.
get.censored = function(i) {
  if (death.times[i] >= 2)
    min(censor.times[i], rexp(1, 0.1))
  else 
    censor.times[i]
}
new.censor.times = sapply(1:N, get.censored)
new.event.times = sapply(1:N, function(i) min(death.times[i], new.censor.times[i]))
new.was.death = death.times < new.censor.times

new.km = gen.km(new.event.times, new.was.death)
plot.km()
lines(new.km, col="blue", do.points=F, xlim=c(0,5))

