from django.shortcuts import render, redirect
from django.http import HttpResponse
from thapp.models import Theory
from thapp.forms import TheoryForm
import os

# Create your views here.
def index(request):
    theories = Theory.objects.using('theory').all()
    context = {'theories': theories}
    return render(request, 'thapp/index.html', context)

def addtheory(request):    
    theory = Theory()    
    if request.method == "POST":
        form = TheoryForm(request.POST, instance=theory)
        if form.is_valid():
            mform = form.save(commit=False)
            mform.save(using='theory')
            os.system('svn ci -m "[thmanager] add new theory to theory.db" theory.db')
            return redirect(index)
    else:
        form = TheoryForm(instance=theory)
    return render(request, 'thapp/addtheory.html', {'form': form, 'theory': theory})
