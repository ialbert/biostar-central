function rev_comp(s) {
    var seq = (new Nt.Seq()).read(s);
    return seq.complement().sequence();
}

function translate(s, offset) {
    offset = typeof offset !== 'undefined' ? offset : 0;
    var seq = (new Nt.Seq()).read(s);
    return seq.translate(offset);
}

function best_match(a, b) {
    var aseq = (new Nt.Seq()).read(a);
    var bseq = (new Nt.Seq()).read(b);
    var aln = aseq.mapSequence(bseq);
    var obj = {
        pos: aln.best().position,
        mask: aln.best().alignmentMask().sequence(),
        aligned: aln.best().alignment().sequence(),
        cover: aln.best().alignmentCover().sequence(),
        obj: aln
    };
    obj.toString = function () {
        return "Best Match: pos=" + this.pos + " mask=" + this.mask;
    }
    return obj;
}


function all_matches(a, b) {
    var aseq = (new Nt.Seq()).read(a);
    var bseq = (new Nt.Seq()).read(b);
    var aln = aseq.mapSequence(bseq);
    return aln;
    var obj = {
        aln: aln
    };
    obj.toString = function () {
        return "All Matches: results=" + aln.results.length + " mask=" + this.mask;
    }
}

function dir(obj) {
    var keys = [];
    for (var key in obj) {
        keys.push(key);
    }
    return keys;
}

function Seq(s) {
    return (new Nt.Seq()).read(s);
}

function efetch(command) {

    var output = 0

    $.ajax(CHARMS_RPC_URL, {
        type: 'POST',
        dataType: 'json',
        async: false,
        data: {command: command},
        success: function (data) {
            if (data.status == 'error') {
                print_error(TERMINAL, data.msg)
            } else {
                try {
                    output = eval(data.result)
                    print_result(TERMINAL, output)
                } catch (e) {
                    print_error(term, e);
                }
            }

        },
        error: function () {
            print_error(TERMINAL, error)
        }
    });

}

OUTPUT_LIMIT = 600
CHARMS_RPC_URL = "/site/charms/rpc/"

function print_result(term, result) {
    if (result !== undefined) {
        out = new String(result);
        if (out.length > OUTPUT_LIMIT) {
            out = out.substring(1, OUTPUT_LIMIT) + "...";
        }
        term.echo(
            "\n" + out + "\n"
        );
    }
}

function print_error(term, e) {
    term.error(new String(e));
}

TERMINAL = 0

$(document).ready(function () {
        $('#charm').terminal(function (command, term) {
            TERMINAL = term;
            if (command !== '') {
                try {
                    var output = window.eval(command);
                    print_result(term, output)
                } catch (e) {
                    term.error(new String(e));
                }
            }
            else {
                term.echo('');
            }
        }, {
            greetings: 'Biostar Charms',
            name: 'charm_terminal',
            height: 500,
            prompt: 'charm> '
        });
    }
);
