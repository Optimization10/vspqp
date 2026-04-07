# Check for Unused Functions

Run the following command to 

```
cppcheck --enable=all --suppress=missingIncludeSystem -I Include -I Source/Include Source Demo 2> /tmp/cppcheck_report.txt
```