module LoggingUtils

using Logging
using LoggingExtras
using Dates
using ProgressMeter

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
    # Date format
    date_format = "yyyy-mm-dd HH:MM:SS"

    # Format function
    function fmt(io, args)
        # [TIME] [LEVEL] [MODULE] Message
        t = Dates.format(now(), date_format)
        l = args.level
        m = args._module
        msg = args.message
        println(io, "[$t] [$l] [$m] $msg")
    end

    # Console Logger
    # We use FormatLogger for consistent formatting.
    console_logger = MinLevelLogger(FormatLogger(fmt, stderr), console_level)

    # File Logger
    # Ensure directory exists
    log_dir = dirname(log_file)
    if !isempty(log_dir) && !isdir(log_dir)
        mkpath(log_dir)
    end

    file_logger = MinLevelLogger(FormatLogger(fmt, log_file; append = true), file_level)

    # Tee Logger
    tee_logger = TeeLogger(console_logger, file_logger)

    # Set global logger
    global_logger(tee_logger)

    @info "Logging initialized. Console level: $console_level, File level: $file_level, Log file: $log_file"
end

end
