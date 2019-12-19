"""Tests for `salt_finder_charts` package."""

from salt_finder_charts.finder_charts import finder_charts


def test_finder_chart_is_a_function():
    """finder_chart is a function."""

    assert callable(finder_charts)
