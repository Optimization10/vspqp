#!/bin/bash

COUNT=3
OUT_FILE="doc/MultipleRunsResults.txt"
OUTPUT_DIR="output"
ARCHIVE_DIR="outputs_archive"

# Append start time
printf "\n\nStart time: %s\n" "$(date)" >> "$OUT_FILE"
printf "Number of runs: %i\n" "$COUNT" >> "$OUT_FILE"

# If a string argument is given, append it first
if [ -n "$1" ]; then
  printf "\n%s\n\n" "$1" >> "$OUT_FILE"
fi

# Suppress normal output but keep error messages
printf "Building project at: %s\n" "$(date)" >> "$OUT_FILE"
make buildPerformance

printf "Running..."
# Run the command in a loop
for i in $(seq 1 $COUNT); do
  # printf "\nRun number %d at %s\n\n" "$i" "$(date)" >> "$OUT_FILE"
  build/vsp >> "$OUT_FILE"

  # If this sleep is not here, then not all runs are saved. It's very
  # inconsistent. What I THINK happens is the following:
  # Bash runs single threaded, and it waits for each command to terminate before
  # it starts the next one. However, it is blind to the file operations
  # initiated in the executable, which the OS might "buffer" meaning they may
  # be completed AFTER the code is done, and therefore after the mv command
  # has been initiated. This is called a race condition.
  sleep 3

  # Archiving step after each run
  if [ -d "$OUTPUT_DIR" ] && [ -n "$(ls -A "$OUTPUT_DIR" 2>/dev/null)" ]; then
    TIMESTAMP=$(date +%m-%d_%H-%M-%S_%Y)
    ARCHIVE_PATH="$ARCHIVE_DIR/Run_$i__$TIMESTAMP"

    # Create archive directory if it doesn't exist
    mkdir -p "$ARCHIVE_PATH"

    # Move files, only report errors if they happen
    if mv "$OUTPUT_DIR"/* "$ARCHIVE_PATH"/; then
      printf "Output files for run %d archived successfully to: %s\n" "$i" "$ARCHIVE_PATH" >> "$OUT_FILE"
    else
      printf "Failed to move files for run %d to %s. Check permissions or directory status.\n" "$i" "$ARCHIVE_PATH" >> "$OUT_FILE"
    fi
  else
    printf "No files to archive for run %d. '%s/' is empty or doesn't exist.\n" "$i" "$OUTPUT_DIR" >> "$OUT_FILE"
  fi
done

# Append end time
printf "\nEnd time: %s\n\n\n" "$(date)" >> "$OUT_FILE"
printf "\nDone!\n"