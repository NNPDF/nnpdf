from django.shortcuts import render, redirect
from django.http import HttpResponse
from thapp.models import Theory
from thapp.forms import TheoryForm
import os

# Create your views here.
def index(request):
    os.system('git pull')
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
            os.system('git add ../../data/theory.db')
            os.system('git commit -m "[apfelcomb] updating theory.db"')
            os.system('git push origin master')
            return redirect(index)
    else:
        form = TheoryForm(instance=theory)
    return render(request, 'thapp/addtheory.html', {'form': form, 'theory': theory})
