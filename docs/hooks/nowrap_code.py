"""Keep every word of inline code on one line.

Browsers break inline code after a hyphen, which splits arguments like `--threads` over two lines. Each
whitespace-separated word of inline code is wrapped in a span that cannot break (see .ausaxs-nowrap in
extra.css), so a long command can still wrap between its words. Code blocks are left alone.
"""

import re

_PRE = re.compile(r"(<pre\b.*?</pre>)", re.S)
_CODE = re.compile(r"<code>([^<]*)</code>")
_WORD = re.compile(r"\S+")


def _wrap(match: re.Match) -> str:
    words = _WORD.sub(lambda w: f'<span class="ausaxs-nowrap">{w.group(0)}</span>', match.group(1))
    return f"<code>{words}</code>"


def on_page_content(html: str, **kwargs) -> str:
    # odd indices are <pre> blocks, which keep their own formatting
    parts = _PRE.split(html)
    return "".join(part if i % 2 else _CODE.sub(_wrap, part) for i, part in enumerate(parts))
