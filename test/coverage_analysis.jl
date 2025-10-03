using Coverage

# process '*.cov' files
coverage = process_folder()
# process '*.info' files, if you collected them
coverage = merge_coverage_counts(
    coverage,
    filter!(
        let prefixes = (joinpath(pwd(), "src", ""), joinpath(pwd(), "deps", ""))
            c -> any(p -> startswith(c.filename, p), prefixes)
        end,
        LCOV.readfolder("test"),
    ),
)
# Get total coverage for all Julia files
covered_lines, total_lines = get_summary(coverage)
