$(document).ready(function () {

    $('.tags > input.search').keydown(function (event) {

        // Prevent submitting form when adding tag by pressing ENTER.
        var ek = event.keyCode || event.which;
        var value = $(this).val().trim();
        if (ek === 13 || ek === 10) {
            event.preventDefault();
            //alert($(this).closest('.tags').html());
            $(this).closest('.tags').dropdown('set selected', value);
            $(this).val('');
            return value
        }
        // Over rides the space
        if (ek === 32) {}
    });

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
