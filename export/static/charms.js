function reverse_complement(s) {
    var seq = (new Nt.Seq()).read(s);
    return seq.complement().sequence();
}

function rc(s) {
    return reverse_complement(s);
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

function __get_default(param, name, default_value) {
    if (typeof param === 'undefined'){
        return default_value
    }
    if (typeof param[name] === 'undefined' ){
        return default_value
    }

    return param[name]
}

function align(query, target, scoring) {
    isLocal = __get_default(scoring, "isLocal", false)
    match = __get_default(scoring, "match",  1)
    misMatch = __get_default(scoring, "misMatch", -1)
    gapOpen = __get_default(scoring, "gapOpen", -1)
    gapExt = __get_default(scoring, "gapExt", -1)
    var aln = bsa_align(isLocal, target, query, [match, misMatch], [gapOpen, gapExt])
    aln['__name'] = 'bsa_align';
    return aln
}
function format(aln){
    var cigar = bsa_cigar2str(aln[2])

    print ('<b>Scores:</b>')
    print('score='+aln[0]+'; pos='+aln[1])
    print('cigar='+ cigar)
    print ('<b>Alignment:</b>')
    var fmt = bsa_cigar2gaps(target, query, aln[1], aln[2])
    print('<code>' + fmt[0] + '</code>')
    print('<code>' + fmt[1] + '</code>')
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
    print("Name=", obj.constructor.name, " attributes:", keys);
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

OUTPUT = $('#charm_output')

function __concat(args) {
    var value = ''
    for (var i = 0; i < args.length; i++) {
        value += ' ' + args[i]
    }
    return value
}

function __append_output(value) {
    OUTPUT.html(OUTPUT.html() + value)
}

function __reset_output() {
    OUTPUT.html("")
}

function clear() {
    __reset_output()
}

function print() {
    __append_output('<div>' + __concat(arguments) + '</div>')
}

function error() {
    __append_output('<div class="error">' + __concat(arguments) + '</div>')
}

function __run(value) {
    try {
        var result = eval(value)
    } catch (e) {
        error(e);
    }
}

function __reset() {
    clear()
    print("<code>Click RUN PROGRAM to run</code>")
}


$(document).ready(function () {

        var editor = ace.edit("editor");
        editor.setTheme("ace/theme/github");
        editor.getSession().setMode("ace/mode/javascript");
        __reset()

        $('#run_charm').click(function () {
            clear();
            var input = editor.getValue();
            __run(input)
        });

        $('.run_me').click(function () {
            var value = $(this).parent().parent().children("pre").text()
            __reset()
            editor.setValue(value)
        });

    }
);
