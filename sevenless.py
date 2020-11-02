#!/usr/bin/env python3

# Simple Sevenless
# Write a program sevenless.py to print out all the numbers from 0 to 99, one on each line, except, do not print any number perfectly divisible by 7.


start = 0
end   = 99
divisor=7
    
print("Printing out numbers from",start,"to",end, " not divisible by",divisor)

for num in range(start,end):
    if num % 7 != 0:
        print(num)
