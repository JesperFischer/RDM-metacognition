

lets fix the x in the inverse logit to 1.

reasons.

we know that the confidence curve shifts with the psychometric function. 

* With x in the numerator this can't happen. (with x = 1 we get the shift of curves we expect).

* The agent doesn't know x so how does he compute confidence using it?

* Fitting the model with x = 1, does nothing to predictive checks for the pilot data (same predictive performance and same estimates for the main thing we are interested in)

* Fitting the model with x = 1, generally just changes meta-bias parameter. In the model with x = 1 we have bias %in% c(-0.1 ; 0.2) in the old its ~ (0.1 ; 0.4). Seems reasonable that not everyone is biased but its a wrong assumption of the model.

* model with x = 1 seems to fit a tiny bit faster

* On the individual subject level the x = 1 model seem to perform better on that we expect (Jesper / Siebe) but worse on data that is perhaps not so expected (Ashley)


