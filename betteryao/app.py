from flask import Flask, render_template, request, session, redirect, url_for, escape
#from OpenSSL import SSL
import subprocess

app = Flask(__name__)

prefsList = []
k = "24" # between 8 and 128
s = "1"
partner_ip = "127.0.0.1"
port = "5880"
mal = "0"
circuit = "gl.4x4x4x4.out.out"

@app.route("/")
def hello():
    return render_template("choose.html")

@app.route("/submit/",methods=['POST'])
def receive_pref():
    if request.method == 'POST':
        print str(request.form)
        uid = int(request.form['id'])
        prefs = request.form['prefs']
        print str(uid)
        print prefs
        print str(len(prefsList))
        prefsList[uid]=prefs
    else:
        pass
    return "receive pref"

@app.route("/create/<numPlayers>/")
def createMatching(numPlayers):
    clear()
    extendMatching(numPlayers)
    return str(numPlayers) + "," + str(len(prefsList))

@app.route("/view/")
def viewMatching():
    return ",".join(str(x) for x in prefsList)

def extendMatching(n):
    prefsList.extend([None for x in xrange(int(n))])

@app.route("/clear/")
def clear():
    prefsList = []
    return str(len(prefsList))

@app.route("/makeInput/<player>/")
def makeInputFile(player):
    pList = ["00000000" if x == None else x for x in prefsList]
    with open("inp."+player+".txt", 'w') as f:
        if player == "gen":
            f.write("".join(str(x)[-3:] for x in pList))
            f.write('\n')
            f.write("".join("000" for x in xrange(len(pList))))
        elif player == "evl":
            f.write("".join("000" for x in xrange(len(pList))))
            f.write('\n')
            f.write("".join(str(x)[-3:] for x in pList))
        else:
            return "error: invalid player"
    return "file written"

@app.route("/runGen/")
def runGen():
    return runMatching("gen")

@app.route("/runEvl/")
def runEvl():
    return runMatching("evl")

def runMatching(player):
    if player in ["gen", "evl"]:
        args = ["./" + player, k, s, circuit, "inp." + player + ".txt", partner_ip, port, mal]
        if mal == "0":
            #return subprocess.check_output(["echo", "Hello " + player + "!"])
            return subprocess.check_output(args)
        elif mal == "1":
            mal_args = ["mpirun","-n",s].extend(args)
            return subprocess.check_output(mal_args)
    else:
        return "error: no valid player called"
    
if __name__ == "__main__":
    app.run(debug=True)
