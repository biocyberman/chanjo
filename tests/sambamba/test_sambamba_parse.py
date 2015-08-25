# -*- coding: utf-8 -*-
from chanjo.sambamba import parse

BASE_HEADERS = ['chrom', 'chromStart', 'chromEnd']
THRESHOLDS = ['percentage10', 'percentage20', 'percentage100']
EXTRA_HEADERS = ['F3', 'F4', 'F5']
SAMBAMBA_HEADERS = ['readCount', 'meanCoverage']

BASE_COLS = ['1', '69089', '70007']
EXTRA_COLS = ['1-69090-70007', '0', '+', 'CCDS30547.1', 'OR4F5']
SAMBAMBA_COLS = ['232', '25.4946']
THRESHOLD_COLS = ['57.9521', '36.0566', '5.55556']
SAMBAMBA_END_COLS = ['ADM992A10', '']


def test_parse_header_no_thresholds():
    header_row = BASE_HEADERS + SAMBAMBA_HEADERS + ['sampleName']
    header = parse.expand_header(header_row)
    assert header['readCount'] == 3
    assert header['thresholds'] == {}
    assert header['extraFields'] == slice(3, 3)
    assert header['sampleName'] == 5


def test_parse_header_single_threshold():
    header_row = (BASE_HEADERS + SAMBAMBA_HEADERS + THRESHOLDS[:1]
                  + ['sampleName'])
    header = parse.expand_header(header_row)
    assert header['readCount'] == 3
    assert header['thresholds'] == {10: 5}
    assert header['sampleName'] == 6


def test_parse_header_multi_thresholds():
    header_row = (BASE_HEADERS + SAMBAMBA_HEADERS + THRESHOLDS
                  + ['sampleName'])
    header = parse.expand_header(header_row)
    assert header['readCount'] == 3
    assert header['thresholds'] == {10: 5, 20: 6, 100: 7}
    assert header['sampleName'] == 8


def test_parse_header_extra_cols():
    header_row = (BASE_HEADERS + EXTRA_HEADERS + SAMBAMBA_HEADERS
                  + ['sampleName'])
    header = parse.expand_header(header_row)
    assert header['readCount'] == 6
    assert header['extraFields'] == slice(3, 6)
    assert header['sampleName'] == 8


def test_parse_row_basic():
    header = {'extraFields': slice(3, 3), 'readCount': 3, 'meanCoverage': 4,
              'thresholds': {}, 'sampleName': 5}
    row = BASE_COLS + SAMBAMBA_COLS + SAMBAMBA_END_COLS
    row_data = parse.expand_row(header, row)
    assert row_data['chrom'] == '1'
    assert row_data['chromStart'] == 69089
    assert row_data['chromEnd'] == 70007
    assert row_data['extraFields'] == []
    assert row_data['readCount'] == 232
    assert row_data['meanCoverage'] == 25.4946
    assert row_data['thresholds'] == {}
    assert row_data['sampleName'] == 'ADM992A10'


def test_parse_row_with_threshold():
    header = {'extraFields': slice(3, 3), 'readCount': 3, 'meanCoverage': 4,
              'thresholds': {10: 5, 20: 6, 100: 7}, 'sampleName': 8}
    row = BASE_COLS + SAMBAMBA_COLS + THRESHOLD_COLS + SAMBAMBA_END_COLS
    row_data = parse.expand_row(header, row)
    assert row_data['thresholds'][10] == 57.9521
    assert row_data['thresholds'][20] == 36.0566
    assert row_data['thresholds'][100] == 5.55556
    assert row_data['sampleName'] == 'ADM992A10'


def test_parse_row_with_extra_fields():
    header = {'extraFields': slice(3, 8), 'readCount': 8, 'meanCoverage': 9,
              'thresholds': {}, 'sampleName': 10}
    row = BASE_COLS + EXTRA_COLS + SAMBAMBA_COLS + SAMBAMBA_END_COLS
    row_data = parse.expand_row(header, row)
    assert row_data['extraFields'] == EXTRA_COLS
    assert row_data['readCount'] == 232
    assert row_data['meanCoverage'] == 25.4946
    assert row_data['sampleName'] == 'ADM992A10'
