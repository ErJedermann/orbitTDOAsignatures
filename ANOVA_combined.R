rec_amounts = c(4,6,8,10,12,14,16,18,20)
m_per_sat = c(1, 3, 5, 7,9, 11, 15, 20, 25, 30)
rec_dia = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
type_sat = c(780, 530)
chunks = 6
false_positive_data = read.csv('./data/experiment1/data_concartenated_both_fp.csv', header = FALSE)
false_negative_data = read.csv('./data/experiment1/data_concartenated_both_fn.csv', header = FALSE)
true_negative_data = read.csv('./data/experiment1/data_concartenated_both_tn.csv', header = FALSE)
true_positive_data = read.csv('./data/experiment1/data_concartenated_both_tp.csv', header = FALSE)
beta = 1/5
response_variable = ((1+beta^2) * true_positive_data) / ((1+beta^2) * true_positive_data + beta^2 * false_negative_data + false_positive_data)
response_variable = array(data = response_variable$V1)
rec_dia2 = rep(rec_dia, times=length(m_per_sat)*length(rec_amounts))  # no direct repetitions, fast cycling (inner axis)
rec_amounts2 = rep(rec_amounts, each=length(rec_dia))  # each value repeats len(rec_dia) times before new cycle begins (middle axis)
rec_amounts2 = rep(rec_amounts2, times=length(m_per_sat))
m_per_sat2 = rep(m_per_sat, each=length(rec_amounts)*length(rec_dia))  # every value repeats len(rec_amounts) * len(rec_dia), no cycles (outer axis)
# multipy the size of each axis due to multiple chunks
rec_dia2 = rep(rec_dia2, times=chunks)
rec_amounts2 = rep(rec_amounts2, times=chunks)
m_per_sat2 = rep(m_per_sat2, times=chunks)
# double the size of each axis due to iridium and starlink
rec_dia2 = rep(rec_dia2, times=2)
rec_amounts2 = rep(rec_amounts2, times=2)
m_per_sat2 = rep(m_per_sat2, times=2)
type_sat = rep(type_sat, each=length(rec_dia)*length(rec_amounts)*length(m_per_sat)*chunks)  # each value repeats len(rec_dia)*len(rec_amounts)*len(m_per_sat)*chunks times before the second value comes
# do the anova
lm_model2 <- lm(response_variable ~ factor(type_sat)*factor(m_per_sat2)*factor(rec_dia2)*factor(rec_amounts2))
anova(lm_model2)

