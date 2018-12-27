# Frequently Asked Questions

## Memory

### Why do I get an `out of memory` error? (answered by Tianjing Zhao)

In Julia, when we run the same code repeatedly, the garbage collection works as follows:


|   Steps    | 1. open| 2. create matrix    | 3. repeat           | 4.repeat            | ...repeat...|
| -----------| -------|---------------------|---------------------|---------------------|-------------|
|   Running  |        |x = rand(10000,10000)|x = rand(10000,10000)|x = rand(10000,10000)|...          |
| Memory(MB) | 89.2   |  866.8              |1628.4               | 1628.4              |1628.4       |

Note that Memory(MB) information are obtained from `Task Manager`(Win10), which shows total physical memory reserved by Julia. You can also check free memory by running `versioninfo(verbose=true) ` or  `Sys.free_memory()/2^20`.

From above table, when we create a matrix, it will definitely use memory, and when we repeat it for the first time (step 3), the memory doubled. But for the consecutive repeating (step 4,5...), the memory doesn't change.

So if you repeat your code for the first time, your will double the memory. **In this situation, restarting kernel will solve the problem.** Here is a **trick** to avoid double memory: change to garbage
(step iii) and collect garbage via `GC.gc()`(step iv) before repeating. See below example:


| Steps      | i. open | ii. create matrix   | iii. set to zero|iv. collect garbage |v. repeat             |
| -----------| --------|---------------------|-----------------|--------------------|----------------------|
| Running    |         |x = rand(10000,10000)| x=0             |GC.gc()             |x = rand(10000,10000) |
| Memory(MB) | 86.9    | 865.1               |865.1            | 98.6               | 861.6                |




## Speed
