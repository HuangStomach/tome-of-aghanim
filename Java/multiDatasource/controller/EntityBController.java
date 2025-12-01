
package com.example.multidatasource.controller;

import com.example.multidatasource.model.b.EntityB;
import com.example.multidatasource.repository.b.EntityBRepository;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.web.bind.annotation.*;

import java.util.List;

@RestController
@RequestMapping("/api/b")
public class EntityBController {
    @Autowired
    private EntityBRepository repository;

    @GetMapping
    public List<EntityB> getAll() {
        return repository.findAll();
    }

    @PostMapping
    public EntityB create(@RequestBody EntityB entity) {
        return repository.save(entity);
    }

    @GetMapping("/{id}")
    public EntityB getById(@PathVariable Long id) {
        return repository.findById(id).orElse(null);
    }

    @PutMapping("/{id}")
    public EntityB update(@PathVariable Long id, @RequestBody EntityB entity) {
        entity.setId(id);
        return repository.save(entity);
    }

    @DeleteMapping("/{id}")
    public void delete(@PathVariable Long id) {
        repository.deleteById(id);
    }
}
