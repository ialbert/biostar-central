function tags_dropdown() {

    $('.tags').dropdown({
        allowAdditions: true,
        // Get form field to add to
        onChange: function (value, text, $selectedItem) {
            // Get form field to add to
            var field = $(this).find("select").data("value");
            var tag_field = $('#{0}'.f(field));
            // Add selected tag to field
            // Set text instead of value
            value = $('<div/>').text(value).html();
            tag_field.val(value);
        }
    });
    $('.tags > input.search').keydown(function (event) {

        // Prevent submitting form when adding tag by pressing ENTER.
        var ek = event.keyCode || event.which;
        var value = $(this).val().trim();

        // Get a list of delimiters
        //var delimiters = $('#field-tags').data('delimiters').split(',');
        // console.log(ek===13);
        if (ek===13) {
            // Escape the text before settings value.
            value = $('<div/>').text(value).html();
            event.preventDefault();
            $(this).closest('.tags').dropdown('set selected', value);
            $(this).val('');
            return value
        }
    })

}

$(document).ready(function () {
    tags_dropdown();
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
