import random
from scipy import stats as scistats


"""
A collection of tests and demos to ensure that new statistical methods and modules are working as intended.

"""

"""
TEST 1: Does the mannwhitneyu function produce the expected output?

list1 = [3,5,1,4,3,5]
list2 = [4,8,6,2,1,9]

null: the two populations are equal
alternative: 'two-sided', the two populations are not equal

expected results:
Value: 13, fail to reject null
"""
print('Mann-U Basic Function Test')
list1 = [3,5,1,4,3,5]
list2 = [4,8,6,2,1,9]
print('list1:', list1)
print('list2:', list2)
print('Expected: Fail to Reject @ p-value of 0.05')
val, pval = scistats.mannwhitneyu(list1, list2, alternative='two-sided')

print('Test Statistic:', val, 'p-value:', pval)
print('Success!')