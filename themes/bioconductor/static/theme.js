


$(document).ready(function () {
    function cancel_answers() {
        $('.answer-text').hide();
        $('.answer-text textarea').val('')
    }

    $('.show-answer').click(function () {
        $(".answer-text").show();
    });

    remove_trigger();

        $('body').keyup(function (event) {
            if (event.keyCode === 27) {
                // Cancel answer text area.
                cancel_answers()
            }
        });

        $('.answer-text .cancel').click(function () {
            cancel_answers()
        });

});
