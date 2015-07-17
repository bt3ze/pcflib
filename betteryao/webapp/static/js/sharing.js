var address1 = "http://127.0.0.1:5000/submit/";
var address2 = "https://www.terner.co";
var words = 4;
//var wordlength = 12;

$(document).ready(function (initfn) {
    console.log("hello world!");
    console.log(initfn);
    sjcl.random.startCollectors();

    $('#submitbtn').click(function () {
        var secretshare1 = generateRandString(words);
        console.log(secretshare1);
        console.log(getSelectedVals());

        var selections = getSelectedVals();
        var secretshare2 = xorPrefs(secretshare1, selections);

        console.log(arrayToString(numsToHexArray(secretshare1)));
        console.log(arrayToString(numsToHexArray(secretshare2)));
        
        sendXORshares(arrayToString(numsToHexArray(secretshare1)),
                        arrayToString(numsToHexArray(secretshare2)));
    });
});

// gets the numerical value
function getSelectedVals() {
    var selects = $('#choose-prefs').find('select.preference option:selected');
    //console.log(selects);
    //_.map(selects, function(a){console.log(a)});
    //console.log(_.map(selects, function(a) { return $(a).val();}));
    return _.map(selects, function (a) {
        return Number($(a).val());
    });
}

/* 
    given a list of selected options and a list of random numbers,
    return a list of the XORs
*/
function xorPrefs(rands, secrets) {
    return _.map(_.zip(rands, secrets), function (o) {
        return o[0] ^ o[1];
    });
}

// hex encodes a list of numbers
function numsToHexArray(ary) {
    return _.map(ary, function (a) {
        if (a < 0) {
            a = a + 0xFFFFFFFF;
        }
        return Number(a).toString(16);
    });
}

function arrayToString(ary) {
    return ary.join("");
}

/*
    address is the place to send the xor-shared inputs
    sharedInput is a string
    cb is a callback function
    cbparams is an array of things to pass to cb
*/
function transmit(address, sharedInput, cb, cbparams) {
    $.ajax({
        url: address,
        data: {
            "prefs": sharedInput,
            "id": 1
        },
        method: "POST",
        success: function (data, textStatus, jqXHR) {
            if (cb) {
                cb.apply(this, cbparams);
            }
        },
        error: function (a, b, c) {
            console.log(a, b, c);
        }
    });
}

// sends both XOR shares off to the server
function sendXORshares(share1, share2) {
    transmit(address1, share1);
    transmit(address2, share2);
}

function generateRandString(numwords) {
    return sjcl.random.randomWords(numwords);
}