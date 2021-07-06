from scripts.ctdna_lengths import tokenize
from itertools import count


def tokens_type_match(types, tokens):
    for type, token in zip(types, tokens):
        assert type == token.type


CHROM = "CHROMOSOME"
NUM = "NUMBER"
TAB = "TAB"
NEWL = "NEWLINE"
COM = "COMMENT"
ANY = "ANY"


class TestTokenize:
    def test_chromosomes(self):
        text = "".join(
            [
                "chr1",
                "chr2",
                "chr3",
                "chr4",
                "chr5",
                "chr6",
                "chr7",
                "chr8",
                "chr9",
                "chr10",
                "chr11",
                "chr12",
                "chr13",
                "chr14",
                "chr15",
                "chr16",
                "chr17",
                "chr18",
                "chr19",
                "chr20",
                "chr21",
                "chr22",
                "chrX",
                "chrY",
            ]
        )
        tokens = list(tokenize(text))
        types = [CHROM] * 24
        for type, token, i in zip(types, tokens, count()):
            assert (
                type == token.type
            ), f"Token value: {token.value}, Tokens successfully processed: {i}"

    def test_tokenize_types_hg38(self):
        text = "chr3\t4000\t5000\nchr3\t6000\t70000"
        tokens = list(tokenize(text))
        types = [CHROM, TAB, NUM, TAB, NUM, NEWL, CHROM, TAB, NUM, TAB, NUM]
        for type, token, i in zip(types, tokens, count()):
            assert (
                type == token.type
            ), f"Token value: {token.value}, Tokens successfully processed: {i}"

    def test_tokenize_types_hg19(self):
        text = "3\t4000\t5000\n3\t6000\t70000"
        tokens = list(tokenize(text))
        types = [NUM, TAB, NUM, TAB, NUM, NEWL, NUM, TAB, NUM, TAB, NUM]
        for type, token, i in zip(types, tokens, count()):
            assert (
                type == token.type
            ), f"Token value: {token.value}, Tokens successfully processed: {i}"
