Description for data in report

To be able to use the input variables as it is done in the (sub-)reports one has to convert the reactives in actual list entries. This is done with the followin line of code
'''
load("..../report/sessionData.RData")
for (n in names(report.env)) {assign(n,report.env[[n]])}
'''

