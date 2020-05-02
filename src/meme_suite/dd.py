# def shortestSubstring(s):
#     set_letter = set(s)
#     if len(set_letter) == len(s):
#         return len(s)
#     ctr = dict()
#     for term in set_letter:
#         ctr[term] = 0
#     start_ptr = 0
#     end_ptr = 0
#     seen_letter = set_letter
#     while seen_letter:
#         letter = s[end_ptr]
#         if not ctr[letter]:
#             ctr[letter] += 1
#             seen_letter.remove(letter)
#         else:
#             ctr[letter] += 1
#         end_ptr += 1
#     min_len = end_ptr - start_ptr
#     print(min_len)
#     while end_ptr != len(s):
#         while True:
#             letter = s[start_ptr]
#             if ctr[letter] == 1:
#                 ctr[letter] -= 1
#                 next_term = letter
#                 start_ptr += 1
#                 break
#             ctr[letter] -= 1
#             start_ptr += 1
#         min_len = min(min_len, end_ptr - start_ptr + 1)
#         while end_ptr != len(s):
#             letter = s[end_ptr]
#             if letter == next_term:
#                 ctr[letter] += 1
#                 end_ptr += 1
#                 break
#             ctr[letter] += 1
#             end_ptr += 1
#     while True:
#         letter = s[start_ptr]
#         if ctr[letter] == 1:
#             ctr[letter] -= 1
#             start_ptr += 1
#             break
#         ctr[letter] -= 1
#         start_ptr += 1
#     min_len = min(min_len, end_ptr - start_ptr + 1)
#     return min_len
# print("ef")
# print(shortestSubstring("bab"))
import re


def countCounterfeit(serialNumber):
    curr_match = dict()
    for term in (10, 20, 50, 100, 200, 500, 1000):
        curr_match[str(term)] = term
    total_value = 0
    for num in serialNumber:
        length = len(num)
        if length > 12 or length < 10:
            continue
        if not re.fullmatch("[A-Z]+", num[:3]):
            continue
        if not re.fullmatch("[A-Z]", num[-1]):
            continue
        if not re.fullmatch('[0-9]+', num[3:7]):
            continue
        year = int(num[3:7])
        if year < 1900 or year > 2019:
            continue
        curr = num[7:-1]
        if curr not in curr_match:
            continue
        total_value += curr_match[curr]

    return total_value


a = ['QDB2012R20B',
'RED190250E',
'RFV201111T',
'TYU20121000E',
'AAA198710B',
'AbC200010E']
print(countCounterfeit(a))