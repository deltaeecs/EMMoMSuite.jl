module LoggingUtils

using Logging
using Dates
using ..LightweightSupport: StreamLogger, TeeLogger, @showprogress

export init_logging, @showprogress

"""
    init_logging(log_file::String; console_level=Logging.Info, file_level=Logging.Debug)

Initialize the global logger with multi-target output:
- Console: Formatted output (default level: Info)
- File: Formatted output (default level: Debug)

# Arguments
- `log_file::String`: Path to the log file.
- `console_level::LogLevel`: Minimum log level for console output.
- `file_level::LogLevel`: Minimum log level for file output.
"""
function init_logging(
    log_file::String;
    console_level::LogLevel = Logging.Info,
    file_level::LogLevel = Logging.Debug,
)
    log_dir = dirname(log_file)
    if !isempty(log_dir) && !isdir(log_dir)
        mkpath(log_dir)
    end

    console_logger = StreamLogger(stderr, console_level)
    file_logger = StreamLogger(open(log_file, "a"), file_level)
    tee_logger = TeeLogger(console_logger, file_logger)

    global_logger(tee_logger)

    @info "Logging initialized. Console level: $console_level, File level: $file_level, Log file: $log_file"
end

end
