#include "shared_api.h"
#include <string.h>
#include <stdarg.h>

#define MAX_LOG_ENTRIES 1E7
#define OUTPUT_DIR "output/logs"

// Global counter for total log entries
static unsigned long log_count = 0;


void Log
(
    const char *filename,
    const char   *format,
    ...
)
{
    // Check if we've reached the maximum number of log entries
    if (log_count >= MAX_LOG_ENTRIES)
    {
        return;  // Silently stop logging
    }
    
    // Create full filepath with output directory
    char filepath[256];
    snprintf(filepath, sizeof(filepath), "%s/%s", OUTPUT_DIR, filename);

    
    // Open file in append mode
    FILE* file = fopen(filepath, "a");
    
    // Write formatted message
    va_list args;
    va_start(args, format);
    vfprintf(file, format, args);
    va_end(args);
    
    // Add newline and close
    fprintf(file, "\n");
    fclose(file);
    
    // Increment counter after successful write
    log_count++;
}