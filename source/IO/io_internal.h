#ifndef IO_INTERNAL_H
#define IO_INTERNAL_H

#include "shared_api.h" // IWYU pragma: keep
#include "io_api.h" // IWYU pragma: keep
#include "mmio.h" // IWYU pragma: keep

VSP_problem *Read_Dot_Ovi_File
(
    const char* filename,
    VSP_options* options
) ;

VSP_problem* VSP_read_MM_file(const char* filename, VSP_options* options) ;

#endif